import sys
import json
import codecs
from utils import encode_num_to_length

call_hex = sys.argv[1]
bytes = bytearray.fromhex(call_hex)

operations = {
  0x01: "Point addition",
  0x02: "Single point multiplication",
  0x03: "Multiexponentiation",
  0x04: "Pairing"
}

pos = 0

operation = bytes[pos:pos+1][0]
pos += 1

def from_be(bs):
  return int(codecs.encode(bs, 'hex'), 16)

print('operation: %s' % operations[operation])

if operation >= 0x01 and operation <= 0x03:
  field_element_length = bytes[pos:pos+1][0]
  pos += 1

  print('field element length: %d' % field_element_length)

  base_field = from_be(bytes[pos:pos+field_element_length])
  print('base field: %d' % base_field)
  pos += field_element_length

  A = from_be(bytes[pos:pos+field_element_length])
  print('A: %d' % A)
  pos += field_element_length

  B = from_be(bytes[pos:pos+field_element_length])
  print('B: %d' % B)
  pos += field_element_length

  scalar_element_length = bytes[pos:pos+1][0]
  pos += 1

  scalar_field_size = from_be(bytes[pos:pos+scalar_element_length])
  print('scalar field size: %d' % scalar_field_size)
  pos += scalar_element_length

  if operation == 0x01:
    x1 = from_be(bytes[pos:pos+field_element_length])
    pos += field_element_length
    y1 = from_be(bytes[pos:pos+field_element_length])
    pos += field_element_length
    print('(x1, y1): (%d, %d)' % (x1, y1))

    x2 = from_be(bytes[pos:pos+field_element_length])
    pos += field_element_length
    y2 = from_be(bytes[pos:pos+field_element_length])
    pos += field_element_length
    print('(x2, y2): (%d, %d)' % (x2, y2))

    F = GF(base_field)
    E = EllipticCurve(F, [F(A), F(B)])
    p1 = E(x1, y1)
    p2 = E(x2, y2)
    p3 = p1 + p2
    print('result: %s' % p3)
    print('result encoded: %s%s' % (encode_num_to_length(p3[0], field_element_length), encode_num_to_length(p3[1], field_element_length)))

  elif operation == 0x02:
    x = from_be(bytes[pos:pos+field_element_length])
    pos += field_element_length
    y = from_be(bytes[pos:pos+field_element_length])
    pos += field_element_length
    print('(x, y): (%d, %d)' % (x, y))

    s = from_be(bytes[pos:pos+scalar_element_length])
    pos += scalar_element_length
    print('s: %d' % s)

    F = GF(base_field)
    E = EllipticCurve(F, [F(A), F(B)])
    p3 = s*E(x, y)
    print('result: %s' % p3)
    print('result encoded: %s%s' % (encode_num_to_length(p3[0], field_element_length), encode_num_to_length(p3[1], field_element_length)))

  elif operation == 0x03:
    pairs = []
    while pos < len(bytes):
      x = from_be(bytes[pos:pos+field_element_length])
      pos += field_element_length
      y = from_be(bytes[pos:pos+field_element_length])
      pos += field_element_length
      s = from_be(bytes[pos:pos+scalar_element_length])
      pos += scalar_element_length
      print('((x, y), s): ((%d, %d), %d)' % (x, y, s))

      pairs.append(((x,y), s))

    F = GF(base_field)
    E = EllipticCurve(F, [F(A), F(B)])

    p = E(0, 1, 0)
    for pair in pairs:
      ((x,y), s) = pair
      p += s*E(x,y)

    print('result: %s' % p)
    print('result encoded: %s%s' % (encode_num_to_length(p[0], field_element_length), encode_num_to_length(p[1], field_element_length)))

else:
  pass
