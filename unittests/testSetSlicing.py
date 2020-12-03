import numc as nc
a = nc.Matrix(3,3)

# testing float cases

a[0][1] = 1.0
print("Expected ", [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], "\n", "Outputted" , a)
print("\n")
a[2][2] = 5
print("Expected ", [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 5.0]], "\n", "Outputted" , a)
print("\n")

b = nc.Matrix(4,4)
b[1:2][2:3] = 5
print("Expected ", [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 5.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]], 
"\n", "Outputted" , b)
b[3:4][3:4] = 8
print("Expected ", [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 5.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 8.0]], 
"\n", "Outputted" , b)
print("\n")
b[2:3][3] = 3
print("Expected ", [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 5.0, 0.0], [0.0, 0.0, 0.0, 3.0], [0.0, 0.0, 0.0, 8.0]], 
"\n", "Outputted" , b)
print("\n")
b[0][1:2] = 4
print("Expected ", [[0.0, 4.0, 0.0, 0.0], [0.0, 0.0, 5.0, 0.0], [0.0, 0.0, 0.0, 3.0], [0.0, 0.0, 0.0, 8.0]], 
"\n", "Outputted" , b)
print("\n")
testMat6 = nc.Matrix(3, 3)
testMat6[0:1, 0:1] = 0.0
testMat6[:, 0] = [1, 1, 1]
expect14 = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
actual14 = testMat6
print("Expected", expect14, "\n", "Outputted" , actual14)

testMat6[0, :] = [2, 2, 2]
expect15 = [[2.0, 2.0, 2.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
actual15 = testMat6
print("Expected", expect15, "\n", "Outputted" , actual15)
