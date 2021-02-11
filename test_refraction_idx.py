import unittest
from refraction_idx import RefractionIdx

class TestRefractionIdx(unittest.TestCase):
    def test_create_refraction_table1(self):
        ref_idx = RefractionIdx()
        table = ref_idx.create_refraction_table(10, 10, 20)
        expected = (10, 10, 20)
        actual = table.shape
        self.assertEqual(expected, actual)

    def test_set_refraction_from_img1(self):
        loader = RefractionIdx()
        loader.set_refraction_from_img("./image/test_pat.tif", 1., 1.2)
        expected_x = 30
        expected_y = 30
        expected_z = 100
        actual_x = loader.num_mesh_x
        actual_y = loader.num_mesh_y
        actual_z = loader.num_mesh_z
        self.assertEqual(expected_x, actual_x)
        self.assertEqual(expected_y, actual_y)
        self.assertEqual(expected_z, actual_z)
        
    def test_set_refraction_from_img(self):
        ref_idx = RefractionIdx()
        ref_idx.set_refraction_from_img()

if __name__ == "__main__":
    unittest.main()