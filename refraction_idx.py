# -*- coding: utf-8 -*-
"""
Created on Mon Jan 4 2021

@author: Kenta Kawaguchi
"""

import cv2
import numpy as np

#### Not Anaconda Library ####
import fdtd

class MeshNum:
    """
    メッシュ数を保持します。

    Attributes
    ----------
    x : int
        X方向のメッシュ数
    y : int
        Y方向のメッシュ数
    z : int
        Z方向のメッシュ数

    """

    def __init__(self, x:int = 0, y:int = 0, z:int = 0):
        self.x = x
        self.y = y
        self.z = z

class RefractionIdx:
    """
    ３次元の屈折率分布を保持します。
    また、各種方法により屈折率分布を生成します。

    Attributes
    -----------
    table : np.ndarray
        屈折率の３次元テーブル
    mesh_num : MeshNum
        メッシュ数。屈折率のテーブルサイズに対応している。

    """
    def __init__(self):
        pass

    def set_refraction_from_img(self, fpath:str, refraction_min:float, refration_max:float):
        """ 
        ２次元画像を読み込んで、GrayScaleに変換し輝度値255を最大屈折率、輝度値0を最小屈折率とする。
        """
        img = cv2.imread(fpath, cv2.IMREAD_GRAYSCALE)
        self.table = self.create_refraction_with_gradation(img, refraction_min, refration_max)

    def create_refraction_with_gradation(self, img: np.ndarray, refraction_min:float, refraction_max:float) -> np.ndarray:
        """ 
        ２次元画像を読み込み、輝度値255を最大屈折率、輝度値0を最小屈折率とする３次元屈折率テーブルを返す。
        """
        self.mesh_num = self.set_mesh_num_from_img(img)
        table = self.create_refraction_table(self.mesh_num)
        for cnt_line, line_img in enumerate(img):
            cross_section = line_img/255*(refraction_max - refraction_min) + refraction_min # xとyの屈折率を同じにする。
            table[:, :, cnt_line] = cross_section
        return table

    def create_refraction_table(self, mesh_num:MeshNum) -> np.ndarray:
        table = np.zeros([mesh_num.x, mesh_num.y, mesh_num.z])
        return table 

    def set_mesh_num_from_img(self, img: np.ndarray) -> MeshNum:
        return MeshNum(img.shape[1], img.shape[1], img.shape[0])

def main():
    loader = RefractionIdx()
    loader.set_refraction_from_img("./image/test_pat.tif", 1., 1.2)
    print(loader.table.shape)
    #print(loader.table)


if __name__ == "__main__":
    main()