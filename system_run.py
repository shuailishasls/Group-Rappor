#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
# @Time    : 19/7/3 20:05
# @Author  : JUN XY
# @Site    : 
# @File    : system_run.py
# @Software: PyCharm

import system
import time


def main():
    start_time = time.time()
    # 初始化系统
    print '系统开始初始化，共用时：%s' % (time.time() - start_time)
    system_create = system.System(200)
    print '系统初始化完成，共用时：%s' % (time.time() - start_time)
    # 设置攻击者人数
    system_create.collect_data(attack_radio=0)
    print '攻击者设置完成，共用时：%s' % (time.time() - start_time)
    # 生成map映射表
    system_create.generate_map()
    print 'map文件生成，共用时：%s' % (time.time() - start_time)
    # 生成少一个组的count文件，与随机count文件
    system_create.minus_one_and_random()
    print 'count文件已生成，共用时：%s' % (time.time() - start_time)
    # 对所有count解码
    system_create.all_count_decode()
    # 将所有解码文件混合
    print '解码文件生成，共用时：%s' % (time.time() - start_time)
    system_create.converge_decode_fit()
    # 计算欧氏距离
    print '混合解码文件完成，共用时：%s' % (time.time() - start_time)
    dis = system_create.euclidean_one_vs_all()
    print dis
    system_create.write_all_data_csv('./csv/test/true.csv')
    print '共用时：%s' % (time.time() - start_time)


if __name__ == '__main__':
    main()
