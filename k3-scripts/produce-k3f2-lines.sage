import time
load("prescribed_roots.sage")

polRing.<x> = PolynomialRing(ZZ)
ans = [polRing(2)]
t = time.time()
c = 0
for i in range(1, 11):
    temp, count = roots_on_unit_circle(2*(x^(2*i)+1), 
                                       filter=no_roots_of_unity,
                                       num_threads=512)
    ans += temp
    c += count
    print len(temp), "polynomials added"
    print c, "nodes enumerated"
    print "time so far: ", time.time() - t, " seconds"

f = open("k3-scripts/k3f2-lines.txt", "wb")
for i in ans:
    f.write(str(i.list()))
    f.write("\n")
f.close()
