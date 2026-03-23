import simple_math as sm

def test_add():
    assert sm.add(3,4) == 7

def test_square_add():
    assert sm.square_add(5,8) == 89

def test_quotrem():
    assert sm.quotrem(11,4)==(2,3)
