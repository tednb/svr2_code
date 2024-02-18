import requests
from bs4 import BeautifulSoup
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import os
import sys
import re
# Function: Get all GSM sample links from a GSE page
def get_gsm_links(gse_url):
    response = requests.get(gse_url)
    soup = BeautifulSoup(response.text, 'html.parser')
    gsm_links = [a['href'] for a in soup.find_all('a', href=True) if 'GSM' in a.text]
    return gsm_links

# Function: Get Characteristics and Title information
def get_characteristics_and_title(gsm_url):
    response = requests.get(gsm_url)
    soup = BeautifulSoup(response.text, 'html.parser')
    characteristics = {}
    
    # Get Characteristics information
    for char in soup.find_all('td', string='Characteristics'):
        char_data = char.find_next_sibling('td').get_text(separator='|').split('|')
        for item in char_data:
            if ':' in item:
                key, value = item.split(':', 1)
                characteristics[key.strip()] = value.strip()
    
    # Get Title information
    title_td = soup.find('td', string='Title')
    if title_td:
        title = title_td.find_next_sibling('td').get_text(strip=True)
        characteristics['Title'] = title

    return characteristics

# Function: Process each GSE series and return the sample count
def process_gse_series(series):
    gse_number = series['GSEnumber']
    gse_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_number}"
    
    # Get all GSM sample links for this GSE series
    gsm_links = get_gsm_links(gse_url)
    
    # Get sample count from the page
    response = requests.get(gse_url)
    soup = BeautifulSoup(response.text, 'html.parser')
    # 查找所有的 <td> 元素
    td_elements = soup.find_all('td')

    # 初始化一个变量来存储找到的样本数
    sample_count = None

    # 遍历所有 <td> 元素
    for td in td_elements:
        if "Samples" in td.get_text():
            # 尝试通过正则表达式获取数字
            match = re.search(r'Samples \((\d+)\)', td.get_text())
            if match:
                sample_count = int(match.group(1))
                break  # 找到后就退出循环

    # 检查是否找到了样本数
    if sample_count is not None:
        print(f"找到的样本数: {sample_count}")
    else:
        print("没有找到包含 'Samples' 的 <td> 元素")

    # Initialize an empty DataFrame to store all sample information for the current GSE series
    gse_samples_data = pd.DataFrame()
    
    # Iterate through each GSM link to get Characteristics information
    for gsm_link in gsm_links:
        full_gsm_url = f"https://www.ncbi.nlm.nih.gov{gsm_link}"
        characteristics = get_characteristics_and_title(full_gsm_url)
        
        # Add GSE number and GSM number to the characteristics dictionary
        characteristics['GSEnumber'] = gse_number
        characteristics['GSMnumber'] = gsm_link.split('=')[1]
        
        # Add the extracted data to the gse_samples_data DataFrame
        gse_samples_data = pd.concat([gse_samples_data, pd.DataFrame([characteristics])], ignore_index=True)
    
    # Save all sample information for the current GSE series to a CSV file
    gse_samples_data_path = os.path.join(output_directory, f'GEO_samples_info_{gse_number}.csv')
    gse_samples_data.to_csv(gse_samples_data_path, index=False)
    print(f'Sample information for GSE series {gse_number} has been saved to file {gse_samples_data_path}.')
    
    return gse_number, sample_count  # Return the GSE number and sample count

# Specify the output directory
output_directory = os.path.expanduser('~/mm285_update')

# If the output directory does not exist, create it
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
    
# List of GEO platform URLs
platform_urls = [
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL32685",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL30650",
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL31950"
]

# Initialize an empty DataFrame to store all results
all_series_data = pd.DataFrame(columns=["GSEnumber", "Title"])

# Iterate through each platform URL
for url in platform_urls:
    response = requests.get(url)
    page = BeautifulSoup(response.text, 'html.parser')

    # Find the outer table containing series information
    outer_tables = page.find_all('table', cellpadding="3")
    
    # Iterate through all outer tables
    for table in outer_tables:
        # Find rows in each outer table
        rows = table.find_all('tr')
        
        # Iterate through all rows to retrieve information
        for row in rows:
            cells = row.find_all('td')
            # Ensure the correct number of cells
            if len(cells) == 2:
                gse_link = cells[0].find('a')
                title = cells[1].get_text(strip=True)
                if gse_link and gse_link.text.startswith('GSE'):
                    gse_number = gse_link.text.strip()

                    # Add the extracted data to the DataFrame
                    new_row = pd.DataFrame([[gse_number, title]], columns=["GSEnumber", "Title"])
                    all_series_data = pd.concat([all_series_data, new_row], ignore_index=True)

geo_series_info_path = os.path.join(output_directory, 'GEO_series_info.csv')
# Check if the file GEO_series_info.csv exists in the current directory
if os.path.exists(geo_series_info_path):
    # Read the existing CSV file content
    existing_data = pd.read_csv(geo_series_info_path)
    
    # Compare the existing data with the newly fetched data
    # Use merge to mark the differences
    merged_data = pd.merge(all_series_data, existing_data, on=["GSEnumber", "Title"], how='outer', indicator=True)

    # Select only the rows that exist in the left DataFrame (i.e., newly fetched data)
    diff_data = merged_data[merged_data['_merge'] == 'left_only']

    # If there are no differences, print a message and exit the program
    if diff_data.empty:
      print("No new data found, the program will not continue running.")
      
    else:
        # If there are differences, only add the differing data to the existing data
        diff_data = diff_data.drop(columns=['_merge'])
        updated_data = pd.concat([existing_data, diff_data], ignore_index=True)
        # Save the updated data to the CSV file
        updated_data.to_csv(geo_series_info_path, index=False)
        print("New data found, file has been updated.")
        # Perform further analysis on the samples under the newly added GSEs
        # Note: Here we only pass the newly added GSE rows to the process_gse_series function
        with ThreadPoolExecutor(max_workers=50) as executor:
            executor.map(process_gse_series, diff_data.to_dict('records'))
      
else:
    # If the CSV file does not exist, save the DataFrame to a CSV
    all_series_data.to_csv(geo_series_info_path, index=False)
    print("CSV file does not exist, new data has been saved, and the program will continue running.")
    # The program continues to run if differences are detected or the file does not exist, the following code will execute
    # Existing DataFrame containing GSEnumber and Title

    all_series_data = pd.read_csv(geo_series_info_path)
    results = []  # Initialize an empty list to store the results

    # Use ThreadPoolExecutor to parallelize the processing of each GSE series
    with ThreadPoolExecutor(max_workers=50) as executor:
        future_results = executor.map(process_gse_series, all_series_data.to_dict('records'))
        results = list(future_results)  # Convert the map object to a list to collect the results

    # Create a new DataFrame to store sample count information
    sample_counts = pd.DataFrame(results, columns=["GSEnumber", "SampleCount"])

    # Merge the sample count information into the all_series_data DataFrame
    all_series_data = all_series_data.merge(sample_counts, on="GSEnumber", how="left")
  
    # Save the DataFrame containing sample count to a CSV file
    geo_series_info_path = os.path.join(output_directory, 'GEO_series_info_with_counts.csv')
    all_series_data.to_csv(geo_series_info_path, index=False)
