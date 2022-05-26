import math
from datetime import datetime
from math import log

def print_curr_time():
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)


def down_sample(in_file_name: str, out_file_name: str, down_to=1000):
    in_file = open(in_file_name, "r")
    count = 0
    for line in in_file:
        count += 1

    in_file.close()
    print(count)

    count_ex_header = count - 1

    decades = log(count_ex_header, 10)

    line_numbers = set([int(math.pow(10, decades/down_to * i)) for i in range(down_to)])

    in_file = open(in_file_name, "r")
    header = in_file.readline()

    out_file = open(out_file_name, "w")
    out_file.write(header)

    count = 0
    for line in in_file:
        if count in line_numbers:
            out_file.write(line)
        count += 1

    in_file.close()
    out_file.close()
