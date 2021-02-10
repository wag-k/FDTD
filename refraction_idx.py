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

    def __init__(self):
        pass

    def set_refraction_from_img(self, fpath:str, refraction_min:float, refration_max:float):
        img = cv2.imread(fpath, cv2.IMREAD_GRAYSCALE)
        self.table = self.create_refraction_with_gradation(img, refraction_min, refration_max)

    def create_refraction_table(self, num_mesh_x: int, num_mesh_y: int, num_mesh_z: int) -> np.ndarray:
        table = np.zeros([num_mesh_x, num_mesh_y, num_mesh_z])
        return table 

    def create_refraction_with_gradation(self, img: np.ndarray, refraction_min:float, refraction_max:float) -> None:
        self.set_mesh_num_from_img(img)
        table = self.create_refraction_table(self.num_mesh_x, self.num_mesh_y, self.num_mesh_z)
        for cnt_line, line_img in enumerate(img):
            cross_section = line_img/255*(refraction_max - refraction_min) + refraction_min
            table[:, :, cnt_line] = cross_section
        return table


    def set_mesh_num_from_img(self, img: np.ndarray) -> None:
        self.num_mesh_x = img.shape[1]
        self.num_mesh_y = img.shape[1]
        self.num_mesh_z = img.shape[0]

def main():
    loader = RefractionIdx()
    loader.set_refraction_from_img("./image/test_pat.tif", 1., 1.2)
    print(loader.table.shape)
    #print(loader.table)


if __name__ == "__main__":
    main()