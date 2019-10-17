# -*- coding:utf-8 -*-
# !/usr/bin/python2
__author__ = 'JUN XY'
__date__ = '19/7/2 14:55'

import encode.rappor as rappor
import random
import numpy as np
import math


class DataGeneration(object):
    def __init__(self, num_strs=100, nonzero=1, expo=10, num_data=1000000):
        self.num_strs = num_strs
        self.nonzero = nonzero
        self.expo = expo
        self.num_data = num_data

    # 生成候选数据集
    def set_of_strings(self):
        result = []
        for i in range(1, self.num_strs + 1):
            result.append('v_' + str(i))

        return result

    # 按指数或线性比例生成数据
    @staticmethod
    def __get_sample_probs(self, decay):
        probs = []
        ind = self.num_strs * self.nonzero
        if decay == "Linear":
            for i in range(1, self.num_strs):
                probs.append(float(i) / sum(range(ind + 1)))
            return probs

        elif decay == "Exponential":
            for i in range(self.num_strs):
                probs.append(float(i) / self.num_strs)
            probs = list(map(lambda x: -x * self.expo, probs[::-1]))
            probs = list(map(lambda x: math.exp(x), probs))
            probs = list(map(lambda x: x / sum(probs), probs))
            return probs[::-1]

        else:
            raise Exception('just can be generate Linear or Exponent data')

    # 生成指定长度指定分布的数据
    def generate_samples(self, decay):
        candidate = self.set_of_strings()
        probs = self.__get_sample_probs(self, decay)

        return np.random.choice(candidate, self.num_data, replace=True, p=probs)

    # 生成攻击者数据集
    def generate_attack(self):
        return


class User(object):
    def __init__(self, user_id, group_number, candidate_num=100):
        self.params = rappor.Params()
        self.id = user_id
        self.upload_data = []
        self.candidate_num = candidate_num
        self.group_number = group_number

    # 指数分布生成数据
    @staticmethod
    def __generation_data(self, num_data, attack=False):
        if attack:
            self.gen_data = ['v_100'] * 10
        else:
            temp = DataGeneration(num_strs=self.candidate_num, num_data=num_data)
            self.gen_data = temp.generate_samples(decay="Exponential").tolist()

    # 数据加扰
    def __perturbation(self, data):
        irr_rand = rappor.SecureIrrRand(self.params)
        enc = rappor.Encoder(self.params, self.group_number, self.id, irr_rand)
        _, prr, _ = enc._internal_encode(data)
        prr_str = rappor.bit_string(prr, self.params.num_bloombits)

        return prr_str

    # 生成上传至服务器的扰动数据
    def upload(self, num_data, attack=False):
        self.__generation_data(self, num_data, attack=attack)
        for data in self.gen_data:
            user_info = [self.id, self.group_number, data, self.__perturbation(data) + '\t']
            self.upload_data.append(user_info)
        return self.upload_data
