import numc as nc
testMat = nc.Matrix(3,3)
expect1 = [0.0,0.0,0.0]
actual1 = testMat[0]
print("Should be", expect1, "it was", actual1)
expect2 = [[0.0,0.0,0.0],[0.0,0.0,0.0]]
actual2 = testMat[0:2]
print("Should be", expect2, "it was", actual2)
expect3 = [[0.0,0.0],[0.0,0.0]]
actual3 = testMat[0:2, 0:2]
print("Should be", expect3, "it was", actual3)
expect4 = [0.0,0.0]
actual4 = testMat[0:2, 0]
print("Should be", expect4, "it was", actual4)
expect5 = [0.0,0.0]
actual5 = testMat[0, 0:2]
print("Should be", expect5, "it was", actual5)
print("Should be", 0.0, "it was", testMat[0,0]) #not sure why this doesnt do 0.0?
testMat2 = nc.Matrix(1,3)
expect6 = 0.0
actual6 = testMat2[0]
print("Should be", expect6, "it was", actual6)
expect7 = [0.0, 0.0]
actual7 = testMat2[0:2]
print("Should be", expect7, "it was", actual7)
testMat3 = nc.Matrix(3,3, 2)
expect8=2.0
actual8 = testMat3[0][1]
print("Should be", expect8, "it was", actual8)
expect9 = 2.0
actual9 = testMat3[0:1,0:1]
print("Should be", expect9, "it was", actual9)
testMat4 = nc.Matrix(3, 1, [1, 2, 3])
expect10 = 1.0
actual10 = testMat4[0]
testMat5 = nc.Matrix(1, 3, [1, 2, 3])
expect11 = 1.0
actual11 = testMat5[0]
expect12 = [2.0, 3.0]
actual12 = testMat4[1:3]
expect13 = [2.0, 3.0]
actual13 = testMat5[1:3]
print("Should be", expect10, "it was", actual10)
print("Should be", expect11, "it was", actual11)
print("Should be", expect12, "it was", actual12)
print("Should be", expect13, "it was", actual13)
testMat6 = nc.Matrix(3, 3)
testMat6[0:1, 0:1] = 0.0
testMat6[:, 0] = [1, 1, 1]
expect14 = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
actual14 = testMat6
testMat6[0, :] = [2, 2, 2]
expect15 = [[2.0, 2.0, 2.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
actual15 = testMat6
testMat6[0:2, 0:2] = [[1, 2], [3, 4]]
expect16 = [[1.0, 2.0, 2.0], [3.0, 4.0, 0.0], [1.0, 0.0, 0.0]]
actual16 = testMat6
print("Should be", expect14, "it was", actual14)
print("Should be", expect15, "it was", actual15)
print("Should be", expect16, "it was", actual16)
