#!/bin/bash

# 下载模型文件
echo "Downloading model..."
curl -L -o /models/model_v1.1.tar.gz 
https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1

# 检查下载是否成功
if [ $? -ne 0 ]; then
    echo "Failed! Check network connection"
    exit 1
fi

# 解压模型文件
echo "unzipping"
tar -xzvf /models/model_v1.1.tar.gz -C /models/

# 检查解压是否成功
if [ $? -ne 0 ]; then
    echo "Failed! Plz redownload it"
    exit 1
fi

echo "Everything finished!"
