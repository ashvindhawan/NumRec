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
testMat6[0:2, 0:2] = [[1, 2], [3, 4]]
expect16 = [[1.0, 2.0, 2.0], [3.0, 4.0, 0.0], [1.0, 0.0, 0.0]]
actual16 = testMat6
print("Expected", expect16, "\n", "Outputted" , actual16)
testMat7 = nc.Matrix(2,2)
testMat7[0:1, 0:1] = 1.0
expect17 = [[1.0,0.0],[0.0,0.0]]
actual17 = testMat7
print("Expected", expect17, "\n", "Outputted" , actual17)
testMat7[1] = [2,2]
expect18 = [[1.0,0.0],[2.0,2.0]]
actual18 = testMat7
print("Expected", expect18, "\n", "Outputted" , actual18)
b = testMat7[1]
b[1] = 3
expect19 = [[1.0,0.0],[2.0,3.0]]
actual19 = testMat7
print("Expected", expect19, "\n", "Outputted" , actual19)
