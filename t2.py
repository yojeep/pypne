import numpy as np
from io import StringIO


st = StringIO("1\t2\t3\t0\t0\n4\t5\t6\n")
array = np.genfromtxt(st,dtype=[('myint','i8'),('myfloat','f8'),('mystring','S5')],delimiter='\t',usecols=(0,1,2))
print(array["myint"])