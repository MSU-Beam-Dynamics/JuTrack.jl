using JuTrack
x = CTPS(0.0, 1, 6, 3) 
px = CTPS(0.0, 2, 6, 3)
y = CTPS(0.0, 3, 6, 3)
py = CTPS(0.0, 4, 6, 3)
z = CTPS(0.0, 5, 6, 3)
delta = CTPS(0.0, 6, 6, 3)

k = 19.6
Q1 = KQUAD(len=0.2, k1=-k, rad=1)
D1 = DRIFT(len=0.2)
Q2 = KQUAD(len=0.2, k1=k, rad=1)
D2 = DRIFT(len=0.2)
B1 = SBEND(len=1.0, angle=0.1, rad=1)

line = [Q1, D1, Q2, D2, B1]
rin = [x, px, y, py, z, delta]

linepass_TPSA!(line, rin, E0=3e9)