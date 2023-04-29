import amplify as amp
import numpy as np

token = 'uQVhv3X4gcNgBOXYIKzau4qqfRdNhZz8'


M = np.array([[1,1,1] , [1,1,0] , [1,0,0]])

def constraints(M):
    a = np.zeros(len(M))
    b = np.zeros(len(M[0]))

    for i in range(len(M)):
        a[i] += np.sum(M[i])
        for k in range(len(M[i])):
            b[k] += M[i,k]

    return a,b


lines,columns = constraints(M)

q = amp.SymbolGenerator(amp.BinaryPoly).array(len(M)*len(M[0]))



print(q)

def Penalty(q, lines, columns, M):

    f = 0

    for i in range(len(lines)):
        for y in range(len(columns)):
            f+= (1 - 2 * lines[i])*q[len(columns)*i + y]

            for z in range(y):
                f+= 2*q[len(columns)*i + y]*q[len(columns)*i + z]


    for i in range(len(columns)):
        for y in range(len(lines)):
            f+= (1 - 2 * columns[i])*q[len(columns)*y + i]

            for z in range(y):
                f+= 2*q[len(columns)*y + i]*q[len(columns)*z + i]


    return f


f = Penalty(q,lines,columns,M)


print(f)


def QUBO(token,f,M):
    client = amp.client.FixstarsClient()
    client.token = token
    client.parameters.timeout = 10000

    model = amp.BinaryQuadraticModel(f)

    solver = amp.Solver(client)
    result = solver.solve(model)

    for solution in result:
        print(f"energy = {solution.energy}\nvalues = {solution.values}")


    for solution in result:
        print(f"q = {q.decode(solution.values)}")

    Q = q.decode(solution.values)
    matrix = M

    for i in range(len(M)):
        for y in range(len(M[0])):
            matrix[i,y] = Q[len(M[0])*i + y]

    return matrix


print(QUBO(token,f,M))
