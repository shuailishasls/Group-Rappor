# -*- coding:utf-8 -*-
# !/usr/bin/python2
__author__ = 'JUN XY'
__date__ = '19/7/2 20:16'

from user import User
from user import DataGeneration
from rpy2.robjects import r
from scipy.spatial import distance
import encode.rappor as rappor
import pandas as pd
import numpy as np
import random
import os

r.source("./decode/decode_dist.R")
os.chdir("../")


def mkdir(path):
    """
    判断路径是否存在，不存在则创建
    :param path: 路径
    """
    is_exists = os.path.exists(path)
    if not is_exists:
        os.makedirs(path)


def iter_files(root_dir):
    """
    遍历多层文件夹下文件
    :param root_dir: 需要变量的文件夹
    :return: 文件夹下的所有文件
    """
    all_file_list = []
    for root, dirs, files in os.walk(root_dir):
        for file_row in files:
            file_name = os.path.join(root, file_row)
            all_file_list.append(file_name)
        for dirname in dirs:
            # 递归调用自身,只改变目录名称
            iter_files(dirname)
    return all_file_list


class System(object):
    def __init__(self, number_user):
        self.users = number_user
        self.__all_data = []
        self.cohort_prr = []
        self.map = []
        self.total_margin = []
        self.root_path = './csv/test'

        mkdir(self.root_path)

        self.params_file = './csv/params.csv'
        self.map_file = os.path.join(self.root_path, 'map.csv')
        self.count_file = os.path.join(self.root_path, 'count')
        self.decode_fit = os.path.join(self.root_path, 'decode_fit')
        self.merge_fit = os.path.join(self.root_path, 'merge_fit.csv')
        self.params = rappor.Params().from_csv(self.params_file)
        self.candidate_string = DataGeneration().set_of_strings()

    def collect_data(self, attack_radio=0):
        """
        采集各个用户所生成的数据
        :param attack_radio: 攻击者的比例
        """
        attack_num = int(self.users * attack_radio)
        value = [int(self.users * 1.0 / self.params.num_cohorts + 0.5)] * self.params.num_cohorts
        all_group = dict(zip(range(0, self.params.num_cohorts), value))

        for num in range(self.users):
            # 初始化攻击者组号
            if num < attack_num:
                group_number = 0
            # 初始化正常用户组号
            else:
                if len(all_group) == 1:
                    group_number = all_group.keys()[0]
                else:
                    group_number = random.choice(all_group.keys())

            all_group[group_number] -= 1
            if all_group[group_number] == 0:
                del all_group[group_number]

            user = User('c' + str(num), group_number)
            user = user.upload(num_data=10, attack=True)
            self.__all_data.extend(user)

        client, cohort, true_value, prr = zip(*self.__all_data)
        self.cohort_prr = zip(cohort, prr)

    @staticmethod
    def write_csv(file_out, data):
        csv_file = pd.DataFrame(data)
        csv_file.to_csv(file_out, header=None, index=False)

    @staticmethod
    def read_csv(file_in, dtype=None, usecols=None, header=None):
        data = pd.read_csv(file_in, dtype=dtype, usecols=usecols, header=header).values.tolist()
        return data

    def write_all_data_csv(self, file_out):
        csv_file = pd.DataFrame(self.__all_data)
        csv_file.to_csv(file_out, header=['client', 'cohort', 'true_value', 'irr'], index=False)

    def __count_margin(self, **data):
        """
        生成边际表数据
        :param data: 需要计算边际表的数据
        :return: 总边际表，没位数量，报告数量
        """
        sums = [[0] * self.params.num_bloombits for _ in range(self.params.num_cohorts)]
        num_reports = [0] * self.params.num_cohorts

        if data:
            data = data['data']
        else:
            data = self.cohort_prr

        for i, row in enumerate(data):
            try:
                (cohort, irr) = row
                irr = irr.strip()
            except ValueError:
                raise RuntimeError('Error parsing row %r' % row)

            if not len(irr) == self.params.num_bloombits:
                raise RuntimeError(
                    "Expected %d bits, got %r" % (self.params.num_bloombits, len(irr)))

            num_reports[int(cohort)] += 1

            for num, c in enumerate(irr):
                if c == '1':
                    sums[cohort][self.params.num_bloombits - num - 1] += 1
                else:
                    if c != '0':
                        raise RuntimeError('Invalid IRR -- digits should be 0 or 1')

        marginal_table = []
        for cohort in range(self.params.num_cohorts):
            marginal_table.append([num_reports[cohort]] + sums[cohort])

        return marginal_table

    def generate_map(self):
        for word in self.candidate_string:
            row = [word]
            for cohort in range(self.params.num_cohorts):
                bloom_bits = rappor.get_bloom_bits(word, cohort, self.params.num_hashes,
                                                   self.params.num_bloombits)
                for bit_to_set in bloom_bits:
                    # bits are indexed from 1.  Add a fixed offset for each cohort.
                    # NOTE: This detail could be omitted from the map file format, and done
                    # in R.
                    row.append(cohort * self.params.num_bloombits + (bit_to_set + 1))
            self.map.append(row)
        self.write_csv(self.map_file, self.map)

    # --------------	 解码	--------------
    def __decode(self, counts_file, file_name):
        """
        解码文件
        :param counts_file: 边际表的地址
        :return: None
        """
        r.Decode_main(self.params_file, self.map_file, counts_file,
                      os.path.join(self.decode_fit, file_name + '_fit.csv'))

    def __combine_fit(self, est_table):
        """
        将解码后生成的fit文件合成一个
        :param est_table: 各fit文件
        :return:
        """
        true = [self.candidate_string]
        for num, file_est in enumerate(est_table):
            estimate = pd.read_csv(file_est, usecols=['string', 'estimate']).values.tolist()
            true.extend([[0] * len(self.candidate_string)])

            for row_est in estimate:
                true[num + 1][int(row_est[0].strip('v_')) - 1] = row_est[-1]

        return [list(row) for row in zip(*true)]

    def all_count_decode(self):
        """
        将所有的count文件解码
        :return: None
        """
        count_file = iter_files(self.count_file)
        mkdir(self.decode_fit)
        for num, file_count in enumerate(count_file):
            file_name = os.path.splitext(os.path.split(file_count)[1])[0]
            self.__decode(file_count, file_name)

    def converge_decode_fit(self):
        """
        将所有的解码文件放入一个表中
        :return: None
        """
        true = self.__combine_fit(iter_files(self.decode_fit))
        pd.DataFrame(true).to_csv(self.merge_fit, header=None, index=None)

    # --------------找攻击者--------------
    def __take_one_cohort(self, num_cohort):
        """
        从总用户上传的数据中取出一个cohort的数据
        :param num_cohort: 需要剔除的组号
        :return: 少一个组的数据
        """
        result = []
        for row in self.cohort_prr:
            if int(row[0]) == num_cohort:
                result.append(row)
        return result

    def __random_minus_count(self):
        """
        总的数据中多次取出与少一个cohort大小相同的数据
        """
        random_1_m = []
        for i in range(self.params.num_cohorts):
            result_one = self.__take_one_cohort(i)
            sample_size = (len(self.cohort_prr) * (self.params.num_cohorts - 1) / (self.params.num_cohorts ** 2))
            result_one_sample = pd.DataFrame(result_one).sample(sample_size, replace=False).values.tolist()
            random_1_m += result_one_sample

        margin_table = self.__count_margin(data=random_1_m)
        self.write_csv(os.path.join(self.count_file, 'random_minus.csv'), data=margin_table)

    def minus_one_and_random(self):
        """
        将count文件提取成m个去掉一个群组的文件与一个随机取出的文件
        :return: None
        """
        mkdir(self.count_file)
        self.total_margin = self.__count_margin()

        for i in range(self.params.num_cohorts):
            result = []
            for j, row in enumerate(self.total_margin):
                if i == j:
                    temp = [0] * (self.params.num_bloombits + 1)
                    result.append(temp)
                else:
                    result.append(row)
            self.write_csv(os.path.join(self.count_file, str(i) + '_minus.csv'), data=result)

        self.__random_minus_count()

    @staticmethod
    def __one_matrix_euclidean(list_2d):
        """
        计算一个矩阵中，每行之间的欧式距离
        :param list_2d: 待计算欧式距离的二维矩阵
        :return: 欧式距离
        """
        dis_result = []
        temp = []
        for i in range(1, len(list_2d) - 1):
            for j in range(i + 1, len(list_2d)):
                temp.extend([i, j, distance.euclidean(list_2d[i], list_2d[j])])
                dis_result.append(temp)
                temp = []
        return dis_result

    def __avg_euclidean(self, euclidean_distance):
        label = 0
        result = [range(1, self.params.num_cohorts + 2), [0] * (self.params.num_cohorts + 1)]
        result = map(lambda x: list(x), zip(*result))

        for row in euclidean_distance:
            for i in range(1, self.params.num_cohorts + 2):
                if i in row:
                    label += 1
                    result[i - 1][-1] += row[-1]
                if label > 1:
                    label = 0
                    break

        result = map(lambda x: x[-1] / (self.params.num_cohorts + 1), result)
        return result

    @staticmethod
    def __find(target, array):
        """
        判断是否存在某个元素
        :param target: 待查找元素
        :param array: 待找寻二维数组
        :return: 如果存在，则返回True，不存在，则返回False
        """
        # 不存在数组则直接返回
        if not array:
            return False
        # 二维数组的行
        row = len(array)
        # 二维数组的列
        col = len(array[0])
        # 二层循环遍历二维数组
        for i in range(row):
            for j in range(col):
                # 如果目标值等于数组中的值，则找到
                if target == array[i][j]:
                    return True
        # 数组遍历结束后仍未找到
        return False

    def euclidean_one_vs_all(self):
        """
        计算merge_fit文件中各列数据间的欧式距离
        :return: 欧氏距离
        """
        data = self.read_csv(self.merge_fit)
        data = [list(row) for row in zip(*data)]
        dis = self.__one_matrix_euclidean(data)
        avg_dis = self.__avg_euclidean(dis)
        avg_dis_std = np.std(avg_dis)
        print avg_dis_std
        return avg_dis_std
