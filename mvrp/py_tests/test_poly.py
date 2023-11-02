import fractions
import numpy as np
from random import randint

px = [randint(1, 9) for _ in range(5)]
px = px + [0 for _ in range(5)]

f_px = np.poly1d(px, r=False, variable=["x"])

z = randint(4,5)
y = f_px(z)
print("z = :", z)
print("y = :", y)

divident = px
divident[-1] -= y
f_divident = np.poly1d(divident, r=False, variable=["x"])

divisor = [1, -z]
f_divisor = np.poly1d(divisor, r=False, variable=["x"])

print(f_divident)
print(f_divisor)

f_qx, rem = np.polydiv(f_divident, f_divisor)
print(f"{f_qx}...{rem}")

r = randint(0, 100)
print(f_qx(r) * f_divisor(r) == f_divident(r))