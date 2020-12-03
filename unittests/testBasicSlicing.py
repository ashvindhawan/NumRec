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
