import unittest
from refraction_idx import RefractionIdx, MeshNum

class TestRefractionIdx(unittest.TestCase):
    def test_create_refraction_table1(self):
        ref_idx = RefractionIdx()
        mesh_num = MeshNum(10, 10, 20)
        table = ref_idx.create_refraction_table(mesh_num)
        expected = (10, 10, 20)
        actual = table.shape
        self.assertEqual(expected, actual)

    def test_set_refraction_from_img1(self):
        loader = RefractionIdx()
        loader.set_refraction_from_img("./image/test_pat.tif", 1., 1.2)
        expected_x = 30
        expected_y = 30
        expected_z = 100
        actual_x = loader.mesh_num.x
        actual_y = loader.mesh_num.y
        actual_z = loader.mesh_num.z
        self.assertEqual(expected_x, actual_x)
        self.assertEqual(expected_y, actual_y)
        self.assertEqual(expected_z, actual_z)
        
    def test_set_refraction_from_img(self):
        ref_idx = RefractionIdx()
        ref_idx.set_refraction_from_img("./image/test_pat.tif", 1., 1.2)
        expected_1 = 1.2
        expected_2 = 1.0
        self.assertEqual(expected_1, ref_idx.table[0, 0, 0])
        self.assertEqual(expected_2, ref_idx.table[29, 29, 99])

if __name__ == "__main__":
    unittest.main()