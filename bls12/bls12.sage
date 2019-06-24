import sys
import json
from utils import encode_num_to_length

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
          ds, R, T, F2, u, E2, RR, F12, w, E12, non_residue, quadratic_non_residue, is_D_type, g1_generator, g2_generator, cofactor_g1, cofactor_g2, A_twist, B_twist, u_to_w = generate_curve(E, b, h, r, x, q, F)
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
          do_pairing(r, x, q, F, ds, R, t, F2, u, E2, RR, F12, w, E12, g1_generator, g2_generator, is_D_type, u_to_w)
          return
      except Exception as e:
        print(e)
        break

def generate_curve(E, b, h, r, x, q, F):
  # common towers from: https://eprint.iacr.org/2012/072.pdf
  R.<T> = PolynomialRing(F)
  non_residue = None
  quadratic_non_residue = None
  F12_equation = None
  u_to_w = None
  if not F(-1).is_square():
    non_residue = -1
    F2.<u> = F.extension(T^2-non_residue,'u')
    for j in range(1,4):
      if not (u+j).is_square():
        quadratic_non_residue = u+j
        F12_equation = (T^6 - j)^2 - non_residue
        u_to_w = T^6 - j
        break
  elif not F(-2).is_square():
    non_residue = -2
    F2.<u> = F.extension(T^2-non_residue,'u')
    if not u.is_square():
      quadratic_non_residue = u
      F12_equation = (T^6)^2 - non_residue
      u_to_w = T^6
    elif not (u+2).is_square():
      quadratic_non_residue = u+2
      F12_equation = (T^6 - 2)^2 - non_residue
      u_to_w = T^6 - 2
  elif not F(-5).is_square():
    non_residue = -5
    F2.<u> = F.extension(T^2-non_residue,'u')
    if not u.is_square():
      quadratic_non_residue = u
      F12_equation = (T^6)^2 - non_residue
      u_to_w = T^6
  if quadratic_non_residue is None:
    raise Exception('can\'t find a quadratic non residue')

  F12.<w> = F.extension(F12_equation)
  print('F12 equation: %s' % F12_equation)
  print('F12 equation is irreducible: %s' % F12_equation.is_irreducible())
  print('non_residue is %s' % non_residue)
  print('quadratic_non_residue is %s' % quadratic_non_residue)
  #print(w^12 % (w^6 + 1))
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

  return ds, R, T, F2, u, E2, RR, F12, w, E12, non_residue, quadratic_non_residue, is_D_type, g1_generator, g2_generator, E.order()/r, E2.order()/r, A_twist, B_twist, u_to_w

def do_pairing(r, x, q, F, ds, R, t, F2, u, E2, RR, F12, w, E12, g1_generator, g2_generator, is_D_type, u_to_w):
  # only works for bls12-377
  #z = miller_loop(E12(81937999373150964239938255573465948239988671502647976594219695644855304257327692006745978603320413799295628339695, 241266749859715473739788878240585681733927191168601896383759122102112907357779751001206799952863815012735208165030), twist(E2(0xea6040e700403170dc5a51b1b140d5532777ee6651cecbe7223ece0799c9de5cf89984bff76fe6b26bfefa6ea16afe*u+0x18480be71c785fec89630a2a3841d01c565f071203e50317ea501f557db6b9b71889f52bb53540274e3e48f7c005196, 0xf8169fd28355189e549da3151a70aa61ef11ac3d591bf12463b01acee304c24279b83f5e52270bd9a1cdd185eb8f93*u+0x690d665d446f7bd960736bcbb2efb4de03ed7274b49a58e458c282f832d204f2cf88886d8c7c2ef094094409fd4ddf), u, w, F2, F12, E2, E12), r, ds, u, w, F2, F12, E2, E12)
  z = miller_loop(E12(g1_generator[0], g1_generator[1]), twist(g2_generator, u, w, F2, F12, E2, E12, is_D_type, u_to_w), r, ds, u, w, F2, F12, E2, E12)
  z2 = miller_loop(E12(g1_generator[0], -g1_generator[1]), twist(g2_generator, u, w, F2, F12, E2, E12, is_D_type, u_to_w), r, ds, u, w, F2, F12, E2, E12)
  z_coeffs = F12.vector_space()(z)
  #print('%s + (%s)*u' % (z_coeffs[0] + z_coeffs[1]*w + z_coeffs[2]*w^2 + z_coeffs[3]*w^3 + z_coeffs[4]*w^4 + z_coeffs[5]*w^5, z_coeffs[6] + z_coeffs[7]*w + z_coeffs[8]*w^2 + z_coeffs[9]*w^3 + z_coeffs[10]*w^4 + z_coeffs[11]*w^5))
  print('pairing result: %s' % ((z*z2)^(int((q^12 - 1)/r))))

def twist(P, u, w, F2, F12, E2, E12, is_D_type, u_to_w):
    x_coeffs = F2.vector_space()(P[0])
    x_multiplier = w^2 if is_D_type else 1/w^2
    y_coeffs = F2.vector_space()(P[1])
    y_multiplier = w^3 if is_D_type else 1/w^3
    return E12((x_coeffs[0] + x_coeffs[1]*u_to_w(w))*x_multiplier, (y_coeffs[0] + y_coeffs[1]*u_to_w(w))*y_multiplier)

def untwist(P, u, w, F2, F12, E2, E12, is_D_type):
    x_coeffs = F12.vector_space()(P[0]/w^2)
    y_coeffs = F12.vector_space()(P[1]/w^3)
    return E2(x_coeffs[0] + x_coeffs[6]*u, y_coeffs[0] + y_coeffs[6]*u)

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


def miller_loop(P, Q, r, ds, u, w, F2, F12, E2, E12):
    #print('in miller loop')
    f = 1
    T = -Q
    L = len(ds)
    for i in range(L-2, 0, -1):
        #print('%d: 1' % i)
        f = f*line_function(T, T, P)
        #print('%d: 2' % i)
        T = 2*T
        #print('%d: 3' % i)
        if ds[i] == -1:
            #print('%d: 4' % i)
            f = f * line_function(T, -Q, P)
            #print('%d: 5' % i)
            T = T - Q
            #print('%d: 6' % i)
        f = f^2

    f = f*line_function(T, T, P)

    return f

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
