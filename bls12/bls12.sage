import sys
import json

search_type = sys.argv[1]
if search_type == 'initial':
  initial_x = int(sys.argv[2], 16)
elif search_type == 'bitsize':
  bitsize = Integer(int(sys.argv[2]))
  x_bitsize = int(bitsize/6)
  initial_x = (2*randrange(0, 2) - 1)*randrange(2**x_bitsize, 2**(x_bitsize+1))
  print('bit size is %d' % bitsize)
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
          ds, R, T, F2, u, E2, RR, TT, F12, w, E12, non_residue, quadratic_non_residue, is_D_type, g1_generator, g2_generator, cofactor_g1, cofactor_g2 = generate_curve(E, b, h, r, x, q, F)
          quadratic_non_residue_coefficients = F2.vector_space()(quadratic_non_residue)
          g2_generator_x_coefficients = F2.vector_space()(g2_generator[0])
          g2_generator_y_coefficients = F2.vector_space()(g2_generator[1])
          g1_test_vectors = generate_scalar_mult_test_vectors_g1(g1_generator, r, 10)
          g2_test_vectors = generate_scalar_mult_test_vectors_g2(g2_generator, F2, r, 10)
          curve_desc = {
            'A': '0',
            'B': num_to_hex(b),
            'x': num_to_hex(x),
            't': num_to_hex(t),
            'q': num_to_hex(q),
            'r': num_to_hex(r),
            'n': num_to_hex(n),
            'non_residue': num_to_hex(non_residue),
            'quadratic_non_residue_0': num_to_hex(quadratic_non_residue_coefficients[0]),
            'quadratic_non_residue_1': num_to_hex(quadratic_non_residue_coefficients[1]),
            'is_D_type': str(is_D_type),
            'g1_x': num_to_hex(g1_generator[0]),
            'g1_y': num_to_hex(g1_generator[1]),
            'g2_x_0': num_to_hex(g2_generator_x_coefficients[0]),
            'g2_x_1': num_to_hex(g2_generator_x_coefficients[1]),
            'g2_y_0': num_to_hex(g2_generator_y_coefficients[0]),
            'g2_y_1': num_to_hex(g2_generator_y_coefficients[1]),
            'cofactor_g1': num_to_hex(cofactor_g1),
            'cofactor_g2': num_to_hex(cofactor_g2),
            'g1_scalar_mult_test_vectors': g1_test_vectors,
            'g2_scalar_mult_test_vectors': g2_test_vectors
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

  return ds, R, T, F2, u, E2, RR, TT, F12, w, E12, non_residue, quadratic_non_residue, is_D_type, g1_generator, g2_generator, E.order()/r, E2.order()/r

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

def generate_scalar_mult_test_vectors_g1(generator, r, amount):
  test_vectors = []
  for i in range(amount):
    b = randrange(0, r)
    a = randrange(0, r)
    p1 = b*generator
    p2 = a*p1

    test_vectors.append({
      'g_x': num_to_hex(p1[0]),
      'g_y': num_to_hex(p1[1]),
      'h_x': num_to_hex(p2[0]),
      'h_y': num_to_hex(p2[1]),
      'a': hex(a),
    });

  return test_vectors

def generate_scalar_mult_test_vectors_g2(generator, F2, r, amount):
  test_vectors = []
  for i in range(amount):
    b = randrange(0, r)
    a = randrange(0, r)
    p1 = b*generator
    p2 = a*p1

    p1_x_coefficients = F2.vector_space()(p1[0])
    p1_y_coefficients = F2.vector_space()(p1[1])
    p2_x_coefficients = F2.vector_space()(p2[0])
    p2_y_coefficients = F2.vector_space()(p2[1])
    test_vectors.append({
      'g_x_0': num_to_hex(p1_x_coefficients[0]),
      'g_x_1': num_to_hex(p1_x_coefficients[1]),
      'g_y_0': num_to_hex(p1_y_coefficients[0]),
      'g_y_1': num_to_hex(p1_y_coefficients[1]),
      'h_x_0': num_to_hex(p2_x_coefficients[0]),
      'h_x_1': num_to_hex(p2_x_coefficients[1]),
      'h_y_0': num_to_hex(p2_y_coefficients[0]),
      'h_y_1': num_to_hex(p2_y_coefficients[1]),
      'a': hex(a),
    });

  return test_vectors

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

def num_to_hex(num):
  if num < 0:
    return '-0x%x' % abs(int(num))
  else:
    return '0x%x' % abs(int(num))

run()
