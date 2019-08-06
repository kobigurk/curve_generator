import sys
import json
from utils import encode_num_to_length

curve = sys.argv[1]
search_type = sys.argv[2]
if search_type == 'initial':
  initial_x = int(sys.argv[3], 16)
elif search_type == 'bitsize':
  bitsize = Integer(int(sys.argv[3]))
  x_bitsize = int(bitsize/4)
  initial_x = (2*randrange(0, 2) - 1)*randrange(2**x_bitsize, 2**(x_bitsize+1))
  print('bit size is %d' % bitsize)
else:
  raise Exception('unknown command')

residue_search_type = sys.argv[4]

out_file = None
if len(sys.argv) > 5:
  out_file = sys.argv[5]

def run():

  for i in range(10**7):
    x = initial_x + i
    if curve == 'bn':
      t = 6*x**2 + 1
      q = 36*x**4 + 36*x**3 + 24*x**2 + 6*x + 1
      r = 36*x**4 + 36*x**3 + 18*x**2 + 6*x + 1
    else:
      R = ZZ['x']
      rx = R.cyclotomic_polynomial(12)
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
        if residue_search_type == 'random':
          b = randrange(-q+1, q)
        E = EllipticCurve(F, [0, b])
        if E.order() == n:
          print('found b: %d' % b)
          print('(x, t, q, r, n): (%d, %d, %d, %d, %d)' % (x, t, q, r, n))
          ds, R, T, F2, u, E2, RR, TT, F12, w, E12, non_residue, quadratic_non_residue, is_D_type, g1_generator, g2_generator, cofactor_g1, cofactor_g2, A_twist, B_twist = generate_curve(E, b, h, r, x, q, F)
          quadratic_non_residue_coefficients = F2.vector_space()(quadratic_non_residue)
          g2_generator_x_coefficients = F2.vector_space()(g2_generator[0])
          g2_generator_y_coefficients = F2.vector_space()(g2_generator[1])

          num_test_vectors = 1
          field_element_length = ceil(log(q)/log(2) * 1/8)
          scalar_element_length = ceil(log(r)/log(2) * 1/8)
          curve_parameters_hex = encode_curve_paramateres(field_element_length, q, 0, b, scalar_element_length, r)
          g1_test_vectors, g1_test_vectors_gs, g1_test_vectors_scalars = generate_scalar_mult_test_vectors_g1(g1_generator, r, num_test_vectors, field_element_length, scalar_element_length, curve_parameters_hex)
          g1_multiexp_test_vector = generate_multiexp_test_vector_g1(E, g1_test_vectors, g1_test_vectors_gs, g1_test_vectors_scalars, field_element_length, scalar_element_length, curve_parameters_hex)
          g2_test_vectors = generate_scalar_mult_test_vectors_g2(g2_generator, F2, r, num_test_vectors)
          A_twist_coefficients = F2.vector_space()(A_twist)
          B_twist_coefficients = F2.vector_space()(B_twist)
          curve_desc = {
            'A': num_to_hex(0),
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
            'g2_scalar_mult_test_vectors': g2_test_vectors,
            'A_twist_0': num_to_hex(A_twist_coefficients[0]),
            'A_twist_1': num_to_hex(A_twist_coefficients[1]),
            'B_twist_0': num_to_hex(B_twist_coefficients[0]),
            'B_twist_1': num_to_hex(B_twist_coefficients[1]),
            'g1_multiexp_test_vector': g1_multiexp_test_vector,
          }
          print('----------------------')
          print(curve_desc)
          if out_file is not None:
            f = open(out_file, 'w')
            f.write(json.dumps(curve_desc))
            f.close()
          return
      except Exception as e:
        print(e)
        break

def generate_curve(E, b, h, r, x, q, F):
  R.<T> = PolynomialRing(F)
  non_residue = None
  quadratic_non_residue = None
  for i in range(10**4):
    num_to_try = None
    if residue_search_type == 'iterate':
      num_to_try = i//2 * (2*(i%2) - 1)
    elif residue_search_type == 'random':
      num_to_try = randrange(-q+1, q)
    else:
      raise Exception('unknown residuce search type: ' + residue_search_type)
    possible_non_residue = F(num_to_try)
    if not possible_non_residue.is_square():
      non_residue = possible_non_residue
    else:
      continue

    F2.<u> = F.extension(T^2-non_residue,'u')

    for j in range(10**2):
      for k in range(10**2):
        num1_to_try = None
        if residue_search_type == 'iterate':
          num1_to_try = j*(j<50) + (j-100)*(j>=50)
        elif residue_search_type == 'random':
          num1_to_try = randrange(-q+1, q)
        else:
          raise Exception('unknown residuce search type: ' + residue_search_type)
        num2_to_try = None
        if residue_search_type == 'iterate':
          num2_to_try = k*(k<50) + (k-100)*(k>=50)
        elif residue_search_type == 'random':
          num2_to_try = randrange(-q+1, q)
        else:
          raise Exception('unknown residuce search type: ' + residue_search_type)

        possible_quadratic_non_residue = F2(num1_to_try + num2_to_try*u)
        if not possible_quadratic_non_residue.is_square():
          quadratic_non_residue = possible_quadratic_non_residue
          break
        else:
          continue
      if quadratic_non_residue is not None:
        break

    if quadratic_non_residue is not None:
      break

  if quadratic_non_residue is None:
    raise Exception('can\'t find a quadratic non residue')

  print('non_residue is %s' % non_residue)
  print('quadratic_non_residue is %s' % quadratic_non_residue)
  ds = Integer(x).digits(2)
  E2 = EllipticCurve(F2, [0,b*quadratic_non_residue])
  is_D_type = False
  A_twist = 0
  if not (E2.order()/r).is_integer():
    is_D_type = true
    E2 = EllipticCurve(F2, [0,b/quadratic_non_residue])
    if not (E2.order()/r).is_integer():
      raise Exception('no twist had appropriate order')
    else:
      B_twist = b/quadratic_non_residue
      print('D type twist')

  else:
    B_twist = b*quadratic_non_residue
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

  return ds, R, T, F2, u, E2, RR, TT, F12, w, E12, non_residue, quadratic_non_residue, is_D_type, g1_generator, g2_generator, E.order()/r, E2.order()/r, A_twist, B_twist

def twist(P, w, E12):
    return E12(P[0]/w^2, P[1]/w^3)

def untwist(P, w, E2):
    return E2((P[0]*w^2)[0], (P[1]*w^3)[0])

def generate_multiexp_test_vector_g1(E, g1_test_vectors, gs, scalars, field_element_length, scalar_element_length, curve_parameters_hex):
  p = E(0, 1, 0)
  str = '03' + curve_parameters_hex
  i = 0
  for v in g1_test_vectors:
    str += encode_num_to_length(int(v['g_x'], 16), field_element_length)
    str += encode_num_to_length(int(v['g_y'], 16), field_element_length)
    str += encode_num_to_length(int(v['a'], 16), scalar_element_length)
    p += scalars[i]*gs[i]

  return {
    'binary': str,
    'expected_x': num_to_hex(p[0]),
    'expected_y': num_to_hex(p[1]),
  }


def generate_scalar_mult_test_vectors_g1(generator, r, amount, field_element_length, scalar_element_length, curve_parameters_hex):
  test_vectors = []
  gs = []
  scalars = []
  for i in range(amount):
    b = randrange(0, r)
    a = randrange(0, r)
    p1 = b*generator
    p2 = a*p1
    p3 = p1 + p2

    gs.append(p1)
    scalars.append(a)
    test_vectors.append({
      'g_x': num_to_hex(p1[0]),
      'g_y': num_to_hex(p1[1]),
      'h_x': num_to_hex(p2[0]),
      'h_y': num_to_hex(p2[1]),
      'a': num_to_hex(a),
      'scalar_mult_binary': '02' + curve_parameters_hex + encode_num_to_length(p1[0], field_element_length) + encode_num_to_length(p1[1], field_element_length) + encode_num_to_length(a, scalar_element_length),
      'addition_binary': '01' + curve_parameters_hex + encode_num_to_length(p1[0], field_element_length) + encode_num_to_length(p1[1], field_element_length) + encode_num_to_length(p2[0], field_element_length) + encode_num_to_length(p2[1], field_element_length),
      'gph_x': num_to_hex(p3[0]),
      'gph_y': num_to_hex(p3[1]),
    });

  return test_vectors, gs, scalars

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


def num_to_hex(num):
  if num < 0:
    return '-0x%x' % abs(int(num))
  else:
    return '0x%x' % abs(int(num))

def encode_curve_paramateres(field_element_length, q, a, b, scalar_element_length, scalar_field_size):
  str = ''
  str += '%x' % field_element_length
  str += encode_num_to_length(q, field_element_length)
  str += encode_num_to_length(a, field_element_length)
  str += encode_num_to_length(b, field_element_length)
  str += '%x' % scalar_element_length
  str += encode_num_to_length(scalar_field_size, scalar_element_length)

  return str

run()
