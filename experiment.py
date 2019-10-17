#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
# @Time	: 19/7/6 10:36
# @Author  : JUN XY
# @Site	:
# @File	: experiment.py
# @Software: PyCharm

import system
import pandas as pd
import numpy as np
import tqdm


def write_data(data, file_out):
    pd.DataFrame(data).to_csv(file_out, index=False, header=None, mode='ab')


def none_10_attack_std():
    """
    测试当攻击者比例从0~0.1时的标准差的变化，其中用户总数为10w个
    """
    data = np.arange(0, 0.11, 0.01)
    for j in tqdm.tqdm(data):
        result = []
        for i in range(50):
            # 初始化系统
            system_create = system.System(100000)
            # 设置攻击者人数
            system_create.collect_data(attack_radio=j)
            # 生成map映射表
            system_create.generate_map()
            # 生成少一个组的count文件，与随机count文件
            system_create.minus_one_and_random()
            # 对所有count解码
            system_create.all_count_decode()
            # 将所有解码文件混合
            system_create.converge_decode_fit()
            # 计算欧氏距离
            dis = system_create.euclidean_one_vs_all()
            result.append(dis)
        result.append('\n')
        write_data(result, './100000_data_attack_0.02_0.1_std.csv')


def attack_01_200_100000_std():
    """
    测试当用户总数为200~10w个，而攻击者比例为0.01时的标准差的变化
    """
    data = range(200, 100200, 200)
    for j in tqdm.tqdm(data):
        result = []
        for i in range(10):
            # 初始化系统
            system_create = system.System(j)
            # 设置攻击者人数
            system_create.collect_data(attack_radio=0.01)
            # 生成map映射表
            system_create.generate_map()
            # 生成少一个组的count文件，与随机count文件
            system_create.minus_one_and_random()
            # 对所有count解码
            system_create.all_count_decode()
            # 将所有解码文件混合
            system_create.converge_decode_fit()
            # 计算欧氏距离
            dis = system_create.euclidean_one_vs_all()
            result.append(dis)
        result.append('\n')
        write_data(result, './200_100000_data_attack_0.01_std.csv')


if __name__ == '__main__':
    attack_01_200_100000_std()
