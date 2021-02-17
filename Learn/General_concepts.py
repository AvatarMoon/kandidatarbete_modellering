import sys
import numpy as np


# Function that yields the product and difference of two
# integers a, b.
# Args:
#     a, b, two integers
#     print_result, bool if result should be printed or not
# Returns:
#     the product and difference of a, b
def mult_and_diff(a, b, print_result=True):
    prod = a * b
    diff = a - b
    if print_result == True:
        # Print float (decimal number) with maximum three decimals
        print("{:.3} * {:.3} = {:.3}".format(a, b, prod))
        print("{:.3} - {:.3} = {:.3}".format(a, b, diff))
    return prod, diff


# Class that holds a string, a list and an integer
class TestClass:
    def __init__(self, a_string, a_list, an_int):
        self.a_string = a_string
        self.a_list = a_list
        self.an_int = an_int


a_list = [1, 2, 3]
a_string = "hello"
a, b = 1, 2
if a == 2 and b != 1:
    print("Case 1")
elif a != 1 or b != 1:
    print("Case 2")
else:
    print("Case 3")

# Print = disp (or fprintf) in Matlab
for char in a_string:
    print(char)
for item in a_list:
    print(item)

# Calling a function 
prod, diff = mult_and_diff(3.3, 4.45)

# Creating a class-object and printing the string
my_obj = TestClass("a_string", [1, 2, 3], 98)
print(my_obj.a_string)

# Some basic numpy
print("Numpy operations")
a = np.array([1.0, 2.0])
b = np.array([1.0, 2.0])
c = a * b   # Multiplies each element (like a .* b)
print("c = ")
print(c)
print("c[0] = {:3f}".format(c[0]))    # Python starts indexing from zero and use bracket-type []

# Operations always acts on entire array
print(np.sqrt(a))
print("Norm of a = {:.3}".format(np.linalg.norm(a, ord=2)))

# Basic matrix operations
a_mat, b_mat = np.array([[1.0, 2], [1, 2]]), np.array([[1.0, 2], [1, 2]])
print(a_mat[0, 0])  # Index within the same bracket
print(a_mat @ b_mat)  # Matrix multiplication
print(a_mat @ a)   # Matrix vector multiplication