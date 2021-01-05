# -*- coding: utf-8 -*-
"""
Created on Mon Jan 4 2021

@author: Kenta Kawaguchi
"""

import cv2
import numpy as np

#### Not Anaconda Library ####
import fdtd

class RefractionIdx:
    def __init__(self, fdtd):
        self.fdtd = fdtd

    def set_refraction_from_img(self, fpath):
        img = cv2.imread(fpath, cv2.IMREAD_GRAYSCALE)
        self.set_mesh_num_from_img(img)
        self.table = self.create_refraction_table(self.num_mesh_x, self.num_mesh_y, self.num_mesh_z)
        # self.set_refraction_with_gradation(img)

    def set_mesh_num_from_img(self, img: np.ndarray) -> None:
        self.num_mesh_x = img.shape[1]
        self.num_mesh_y = img.shape[1]
        self.num_mesh_z = img.shape[0]
        
    def set_refraction_with_gradation(self, img: np.ndarray) -> None:
        self.fdtd.diff_table = np.zeros([self.fdtd.n_meshx, self.fdtd.n_meshy, self.fdtd.n_meshz])

    def create_refraction_table(self, num_mesh_x: int, num_mesh_y: int, num_mesh_z: int) -> np.ndarray:
        table = np.zeros([num_mesh_x, num_mesh_y, num_mesh_z])
        return table


def main():
    fdtd_ = fdtd.FDTD()
    loader = RefractionIdx(fdtd_)
    loader.set_refraction_from_img("./image/test_pat.tif")
    print(loader.table)


if __name__ == "__main__":
    main()