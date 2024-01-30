import unittest
import tpsa

class TestTPSA(unittest.TestCase):
    da_dim = 3
    da_order = 4
    n_vec = 4000
    eps = 1e-14
    path = "../tpsa/test/"

    def test_00_DAIni(self):  
        init = tpsa.da_init(self.da_order, self.da_dim, self.n_vec)
        self.assertEqual(init, 0)
        
    def test_01_Exp(self):
        da = tpsa.base()
        x = 1 + da[0] + 2*da[1] + 5*da[2]
        y = tpsa.exp(x)
        self.assertTrue(tpsa.compare_da_with_file(self.path+"exp_da.txt", y, self.eps))
        
    def test_02_SubstNum(self):
        da = tpsa.base()
        x = 1 + da[0] + 2*da[1] + 5*da[2]
        y = tpsa.exp(x)
        z = tpsa.assign()
        tpsa.da_substitute_const(y, 0, 1, z);
        self.assertTrue(tpsa.compare_da_with_file(self.path+"substitute_number.txt", z, self.eps))
        
    def test_03_SubstDA(self):
        da = tpsa.base()
        x = 1 + da[0] + 2*da[1] + 5*da[2]
        y = tpsa.exp(x)
        z = tpsa.assign()
        tpsa.da_substitute(y, 0, x, z);
        self.assertTrue(tpsa.compare_da_with_file(self.path+"substitute_da_vector.txt", z, self.eps))
        
    def test_04_SubstDAs(self):
        da = tpsa.base()
        x = 1 + da[0] + 2*da[1] + 5*da[2]
        y = tpsa.exp(x)
        z = tpsa.assign()
        idx = [0,1]
        lv = tpsa.assign(2)
        lv[0] = tpsa.sin(x)
        lv[1] = tpsa.cos(x)
        tpsa.da_substitute(y, idx, lv, z);
        self.assertTrue(tpsa.compare_da_with_file(self.path+"substitute_multiple_da_vectors.txt", z, self.eps))
        
    def test_05_SubstBunch(self):
        da = tpsa.base()
        x = 1 + da[0] + 2*da[1] + 5*da[2]
        y = tpsa.exp(x)
        z = tpsa.assign()
        idx = [0,1]
        lv = tpsa.assign(2)
        lv[0] = tpsa.sin(x)
        lv[1] = tpsa.cos(x)
        lx = tpsa.assign(3)
        lx[0] = x
        lx[1] = y
        lx[2] = tpsa.sinh(x)
        ly = tpsa.assign(3)
        tpsa.da_substitute(lx, idx, lv, ly);
        self.assertTrue(tpsa.compare_da_with_file(self.path+"bunch_substitution_0.txt", ly[0], self.eps))
        self.assertTrue(tpsa.compare_da_with_file(self.path+"bunch_substitution_1.txt", ly[1], self.eps))
        self.assertTrue(tpsa.compare_da_with_file(self.path+"bunch_substitution_2.txt", ly[2], self.eps))
        
    def test_06_Composition(self):
        da = tpsa.base()
        x = 1 + da[0] + 2*da[1] + 5*da[2]
        y = tpsa.exp(x)
        lx = tpsa.assign(3)
        lx[0] = x
        lx[1] = y
        lx[2] = tpsa.sinh(x)
        lu = tpsa.assign(3)
        lu[0] = tpsa.sin(x)
        lu[1] = tpsa.cos(x)
        lu[2] = tpsa.tan(x)
        ly = tpsa.assign(3)
        tpsa.da_composition(lx, lu, ly)
        self.assertTrue(tpsa.compare_da_with_file(self.path+"da_composition_0.txt", ly[0], self.eps))
        self.assertTrue(tpsa.compare_da_with_file(self.path+"da_composition_1.txt", ly[1], self.eps))
        self.assertTrue(tpsa.compare_da_with_file(self.path+"da_composition_2.txt", ly[2], self.eps))
        
    def test_07_CDfuncs(self):
        da = tpsa.base()
        x1 = da[0] + 2*da[1] + 3*da[2]
        x2 = tpsa.sin(x1)
        x1 = tpsa.cos(x1)
        x3 = 0.5*da[0] + 4*da[1] + 2.7*da[2]
        x4 = tpsa.sin(x3)
        x3 = tpsa.cos(x3)
        y1 = tpsa.complex(x1, x2)
        y2 = tpsa.complex(x3, x4)
        r = y1+y2
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"cd_calculation_0.txt", r, self.eps))
        r = y1-y2
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"cd_calculation_1.txt", r, self.eps))
        r = y1*y2
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"cd_calculation_2.txt", r, self.eps))
        r = y1/y2
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"cd_calculation_3.txt", r, self.eps))
        
    def test_08_DAComposCD(self):
        da = tpsa.base()
        x1 = da[0] + 2*da[1] + 3*da[2]
        x2 = tpsa.sin(x1)
        x1 = tpsa.cos(x1)
        x3 = 0.5*da[0] + 4*da[1] + 2.7*da[2]
        x4 = tpsa.sin(x3)
        x3 = tpsa.cos(x3)
        y1 = tpsa.complex(x1, x2)
        y2 = tpsa.complex(x3, x4)
        mmap = tpsa.assign(2)
        mmap[0] = x1
        mmap[1] = x2
        cnmap = tpsa.assign_cd(3)
        cnmap[0]= y1
        cnmap[1] = y2
        cnmap[2] = y1*y2
        comap = tpsa.assign_cd(2)
        tpsa.cd_composition(mmap, cnmap, comap)
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"da_composition_cd_0.txt", comap[0], self.eps))
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"da_composition_cd_1.txt", comap[1], self.eps))
        
    def test_09_CDComposCD(self):
        da = tpsa.base()
        x1 = da[0] + 2*da[1] + 3*da[2]
        x2 = tpsa.sin(x1)
        x1 = tpsa.cos(x1)
        x3 = 0.5*da[0] + 4*da[1] + 2.7*da[2]
        x4 = tpsa.sin(x3)
        x3 = tpsa.cos(x3)
        y1 = tpsa.complex(x1, x2)
        y2 = tpsa.complex(x3, x4)
        cmmap = tpsa.assign_cd(2)
        cmmap[0] = tpsa.complex(x1,tpsa.exp(x1))
        cmmap[1] = tpsa.complex(x2,tpsa.exp(x2))
        cnmap = tpsa.assign_cd(3)
        cnmap[0]= y1
        cnmap[1] = y2
        cnmap[2] = y1*y2
        comap = tpsa.assign_cd(2)
        tpsa.cd_composition(cmmap, cnmap, comap)
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"cd_composition_cd_0.txt", comap[0], self.eps))
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"cd_composition_cd_1.txt", comap[1], self.eps))
        
    def test_09_CDComposDA(self):
        da = tpsa.base()
        x1 = da[0] + 2*da[1] + 3*da[2]
        x2 = tpsa.sin(x1)
        x1 = tpsa.cos(x1)
        x3 = 0.5*da[0] + 4*da[1] + 2.7*da[2]
        x4 = tpsa.sin(x3)
        x3 = tpsa.cos(x3)
        y1 = tpsa.complex(x1, x2)
        y2 = tpsa.complex(x3, x4)
        cmmap = tpsa.assign_cd(2)
        cmmap[0] = tpsa.complex(x1,tpsa.exp(x1))
        cmmap[1] = tpsa.complex(x2,tpsa.exp(x2))
        mmap = tpsa.assign(3)
        mmap[0] = x1
        mmap[1] = x2
        mmap[2] = x1+0.33*x2
        comap = tpsa.assign_cd(2)
        tpsa.cd_composition(cmmap, mmap, comap)
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"cd_composition_da_0.txt", comap[0], self.eps))
        self.assertTrue(tpsa.compare_cd_with_file(self.path+"cd_composition_da_1.txt", comap[1], self.eps))
    
if __name__ == '__main__':
    unittest.main()