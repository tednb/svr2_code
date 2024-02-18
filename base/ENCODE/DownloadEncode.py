import requests
import os
import concurrent.futures

# ENCODE的API的基本URL
base_url = "https://www.encodeproject.org"

# 你的实验列表
experiment_ids = ['ENCSR735PMJ','ENCSR019HJV','ENCSR265OMO','ENCSR227RPE','ENCSR842QTB','ENCSR906HLA','ENCSR695BII','ENCSR324F','ENCSR027ICI','ENCSR215QTA','ENCSR476OEL','ENCSR013VIR','ENCSR430SAY','ENCSR050EXR','ENCSR037RNN','ENCSR217TMK','ENCSR425NDU','ENCSR129SBE','ENCSR535YCH','ENCSR871FJZ','ENCSR953WFU'
]  # 你可以添加更多的实验ID

def download_file(file_info, experiment_id):
    # 检查文件类型、状态和版本
    # if (file_info.get("file_type") == "tsv" and 
    #     file_info.get("output_type") == "gene quantifications" and 
    #     file_info.get("status") == "released" and 
    #     file_info.get("genome_annotation") == "M21"): 
      if (file_info.get("file_type") == "bed bed9+" and 
        file_info.get("output_type") == "methylation state at CpG" and 
        file_info.get("status") == "released"): 
        # 获取文件的下载URL
        download_url = file_info["href"]

        # 创建以实验ID命名的文件夹
        directory = f"data//MouseDNAmAtlas/ENCODE/Embryonic/WGBS/E_16/{experiment_id}"
        os.makedirs(directory, exist_ok=True)

        # 获取文件名
        file_name = os.path.join(directory, download_url.split("/")[-1])

        # 下载文件
        response = requests.get(f"{base_url}{download_url}", stream=True)

        # 将文件写入磁盘
        with open(file_name, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        print(f"File {file_info['accession']} downloaded")

with concurrent.futures.ThreadPoolExecutor() as executor:
    for experiment_id in experiment_ids:
        # 获取实验信息
        response = requests.get(f"{base_url}/experiments/{experiment_id}/?format=json")
        experiment_info = response.json()

        # 检查实验是否存在
        if experiment_info.get("status") != "released":
            print(f"Experiment {experiment_id} not found or not available")
            continue

        # 遍历实验中的文件
        for file_info in experiment_info.get("files", []):
            # 使用线程池下载文件
            executor.submit(download_file, file_info, experiment_id)
