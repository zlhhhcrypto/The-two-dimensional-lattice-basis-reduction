import random
import sys
from decimal import Decimal, getcontext

getcontext().prec = 2000000

sys.set_int_max_str_digits(2000000)

def generate_random_decimal(num_digits):
    if num_digits < 1:
        raise ValueError("0")
    first_digit = random.randint(1, 9)
    remaining_digits = [random.randint(0, 9) for _ in range(num_digits - 1)]
    random_number = int(str(first_digit) + ''.join(map(str, remaining_digits)))

    return random_number

num_digits = 200000
ay = generate_random_decimal(num_digits)
by = generate_random_decimal(num_digits)

if ay < by:
    ay, by = by, ay

def fraction_to_continued_fraction(numerator, denominator):
    if denominator == 0:
        raise ValueError("Denominator cannot be zero")

    continued_fraction = []
    while denominator != 0:
        quotient = numerator // denominator
        continued_fraction.append(quotient)
        numerator, denominator = denominator, numerator % denominator

    return continued_fraction

continued_fraction_ay_by = fraction_to_continued_fraction(ay, by)
continued_fraction_ax_bx = fraction_to_continued_fraction(ay, by)  # 初始时让两个连分数展开一样

#half_length = len(continued_fraction_ay_by)-1
#half_length = len(continued_fraction_ay_by)*7 // 8
#half_length = len(continued_fraction_ay_by)*6 // 8

#half_length = len(continued_fraction_ay_by)*29 // 40
#half_length = len(continued_fraction_ay_by)*7 // 10

#half_length = len(continued_fraction_ay_by)*5 // 8
#half_length = len(continued_fraction_ay_by)*4 // 8
#half_length = len(continued_fraction_ay_by)*3 // 8
#half_length = len(continued_fraction_ay_by)*2 // 8
#half_length = len(continued_fraction_ay_by)*1 // 8
#half_length = len(continued_fraction_ay_by)*1 // 20
#half_length = len(continued_fraction_ay_by)*1 // 40
#half_length = len(continued_fraction_ay_by)*1 // 200
#half_length = len(continued_fraction_ay_by)*1 // 400
#half_length = len(continued_fraction_ay_by)*1 // 2000
#half_length = len(continued_fraction_ay_by)*1 // 4000
#half_length = len(continued_fraction_ay_by)*1 // 20000
#half_length = len(continued_fraction_ay_by)*1 // len(continued_fraction_ay_by)
#half_length = len(continued_fraction_ay_by)*29 // 40
#half_length = len(continued_fraction_ay_by)
for i in range(half_length):
    continued_fraction_ax_bx[i] = continued_fraction_ay_by[i]
for i in range(half_length, len(continued_fraction_ay_by)):
    if continued_fraction_ay_by[i] > 0:
    	continued_fraction_ax_bx[i] = continued_fraction_ay_by[i] + random.randint(-1, 0)
    else:
    	continued_fraction_ax_bx[i] = continued_fraction_ay_by[i]

def continued_fraction_to_fraction(continued_fraction):
    if not continued_fraction:
        raise ValueError("Continued fraction cannot be empty")

    numerator, denominator = 1, 0
    for value in reversed(continued_fraction):
        numerator, denominator = value * numerator + denominator, numerator

    return numerator, denominator

ax, bx = continued_fraction_to_fraction(continued_fraction_ax_bx)
if ax < 0:
    ax = -ax
    bx = -bx
with open('Input.txt', 'w') as file:
    file.write(f"{ay}\n")
    file.write(f"{ax}\n")
    file.write(f"{by}\n")
    file.write(f"{bx}")