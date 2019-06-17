import sys
import json

search_type = sys.argv[1]
if search_type == 'initial':
  initial_x = int(sys.argv[2], 16)
elif search_type == 'bitsize':
  bitsize = Integer(int(sys.argv[2]))
  x_bitsize = int(bitsize/6)
  initial_x = randrange(2**x_bitsize, 2**(x_bitsize+1))
else:
  raise Exception('unknown command')

out_file = None
if len(sys.argv) > 3:
  out_file = sys.argv[3]

def run():
  R = ZZ['x']
  rx = R.cyclotomic_polynomial(12)

  for i in range(10**7):
    x = initial_x + i
    t = x + 1
    q = 1/3 * (x - 1)**2*(x**4 - x**2 + 1) + x
    r = rx(x)
    if q.is_integer():
      q = int(q)
    else:
      continue
    n = q + 1 - t
    if not is_prime(q):
      continue
    h = int(n/r)
    #print('prime q: %d' % q)

    for b in range(1, 10**7):
      try:
        F = GF(q)
        E = EllipticCurve(F, [0, b])
        if E.order() == n:
          print('found b: %d' % b)
          print('(x, t, q, r, n): (%d, %d, %d, %d, %d)' % (x, t, q, r, n))
          ds, R, T, F2, u, E2, RR, TT, F12, w, E12, non_residue, quadratic_non_residue, is_D_type, g1_generator, g2_generator = generate_curve(E, b, h, r, x, q, F)
          quadratic_non_residue_coefficients = F2.vector_space()(quadratic_non_residue)
          g2_generator_x_coefficients = F2.vector_space()(g2_generator[0])
          g2_generator_y_coefficients = F2.vector_space()(g2_generator[1])
          curve_desc = {
            'x': hex(int(x)),
            't': hex(int(t)),
            'q': hex(int(q)),
            'r': hex(int(r)),
            'n': hex(int(n)),
            'non_residue': hex(int(non_residue)),
            'quadratic_non_residue_0': hex(int(quadratic_non_residue_coefficients[0])),
            'quadratic_non_residue_1': hex(int(quadratic_non_residue_coefficients[1])),
            'is_D_type': str(is_D_type),
            'g1_x': hex(int(g1_generator[0])),
            'g1_y': hex(int(g1_generator[1])),
            'g2_x_0': hex(int(g2_generator_x_coefficients[0])),
            'g2_x_1': hex(int(g2_generator_x_coefficients[1])),
            'g2_y_0': hex(int(g2_generator_y_coefficients[0])),
            'g2_y_1': hex(int(g2_generator_y_coefficients[1]))
          }
          print('----------------------')
          print(curve_desc)
          if out_file is not None:
            f = open(out_file, 'w')
            f.write(json.dumps(curve_desc))
            f.close()
          #do_pairing(r, x, q, F, ds, R, t, F2, u, E2, RR, tt, F12, w, E12)
          return
      except Exception as e:
        print(e)
        break

def generate_curve(E, b, h, r, x, q, F):
  # common towers from: https://eprint.iacr.org/2012/072.pdf
  R.<T> = PolynomialRing(F)
  non_residue = None
  quadratic_non_residue = None
  if not F(-1).is_square():
    non_residue = -1
    F2.<u> = F.extension(T^2-non_residue,'u')
    for j in range(1,4):
      if not (u+j).is_square():
        quadratic_non_residue = u+j
        break
  elif not F(-2).is_square():
    non_residue = -2
    F2.<u> = F.extension(T^2-non_residue,'u')
    if not u.is_square():
      quadratic_non_residue = u
    elif not (u+2).is_square():
      quadratic_non_residue = u+2
  elif not F(-5).is_square():
    non_residue = -5
    F2.<u> = F.extension(T^2-non_residue,'u')
    if not u.is_square():
      quadratic_non_residue = u
  if quadratic_non_residue is None:
    raise Exception('can\'t find a quadratic non residue')

  print('non_residue is %s' % non_residue)
  print('quadratic_non_residue is %s' % quadratic_non_residue)
  ds = Integer(x).digits(2)
  E2 = EllipticCurve(F2, [0,b*quadratic_non_residue])
  is_D_type = False
  if not (E2.order()/r).is_integer():
    is_D_type = true
    E2 = EllipticCurve(F2, [0,b/quadratic_non_residue])
    if not (E2.order()/r).is_integer():
      raise Exception('no twist had appropriate order')
    else:
      print('D type twist')

  else:
    print('M type twist')

  RR.<TT> = PolynomialRing(F2)
  F12.<w> = F2.extension(TT^6 - quadratic_non_residue)
  E12 = EllipticCurve(F12, [0,b])

  g1_generator = None
  for j in range(10**3):
    attempted_x = F(j)
    y2 = attempted_x**3 + F(b)
    if not y2.is_square():
      continue
    y = y2.sqrt()
    p = h*E(attempted_x, y)
    if not p.is_zero() and (r*p).is_zero():
      g1_generator = p
      print('found generator for G1: %s' % p)
      break

  g2_generator = None
  h2 = int(E2.order()/r)
  for j in range(10**3):
    attempted_x = F2(j)
    if is_D_type:
      b2 = b/quadratic_non_residue
    else:
      b2 = b*quadratic_non_residue

    y2 = attempted_x**3 + F2(b2)
    if not y2.is_square():
      continue
    y = y2.sqrt()
    p = h2*E2(attempted_x, y)
    if not p.is_zero() and (r*p).is_zero():
      g2_generator = p
      print('found generator for G2: %s' % p)
      break

  return ds, R, T, F2, u, E2, RR, TT, F12, w, E12, non_residue, quadratic_non_residue, is_D_type, g1_generator, g2_generator

def do_pairing(r, x, q, F, ds, R, t, F2, u, E2, RR, tt, F12, w, E12):
  # only works for bls12-381
  #z = miller_loop(E12(3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507, 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569), twist(E2(3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758*u+352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160, 927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582*u+1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905), w, E12), r, ds, w, E2, E12)
  z = miller_loop(E12(3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507, 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569), twist(E2(3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758*u+352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160, 927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582*u+1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905), w, E12), r, ds, w, E2, E12)
  print(z)
  print(z^(int((q^12 - 1)/r)))

def twist(P, w, E12):
    return E12(P[0]/w^2, P[1]/w^3)

def untwist(P, w, E2):
    return E2((P[0]*w^2)[0], (P[1]*w^3)[0])

def line_function(A, B, P):
    if A==B:
        l = (3*A[0]^2)/(2*A[1])
    elif A == -B:
        return P[0]-A[0]
    else:
        l = (B[1]-A[1])/(B[0]-A[0])
    return l*(P[0]-A[0]) + A[1] - P[1]


def miller_loop(P, Q, r, ds, w, E2, E12):
    f = 1
    T = -Q
    L = len(ds)
    for i in range(L-2, 0, -1):
        print('%d: 1' % i)
        f = f*line_function(T, T, P)
        print('%d: 2' % i)
        T = twist(2*untwist(T, w, E2), w, E12)
        print('%d: 3' % i)
        if ds[i] == -1:
            print('%d: 4' % i)
            f = f * line_function(T, -Q, P)
            print('%d: 5' % i)
            T = twist(untwist(T, w, E2) - untwist(Q, w, E2), w, E12)
            print('%d: 6' % i)
        f = f^2

    f = f*line_function(T, T, P)

    return f

run()
