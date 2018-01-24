import pytest

POTENTIAL = '''
#F 3 3
#G 2 3 3
#I 0 0 0
#E

1.7785714285714285e+00 9.0000000000000000e+00 14
4.6483646502270971e-02 7.3000328284430382e-01 6
{} {} 6

-2.1847554330621324e+00 0.0000000000000000e+00
8.4255614529302880e-01
2.0271471509245159e-01
2.2155522312087100e-02
-6.8303638148139376e-02
-4.9043278969706922e-02
-1.0946697911942214e-02
-1.0298525723152156e-02
-6.5685589306463460e-03
-2.0283626156443295e-03
1.2197714420628298e-02
2.2148801387662245e-02
1.4461939319528331e-02
9.9020972820913523e-03
0.0000000000000000e+00

-9.7346191762795655e+00 -8.3400770422903880e-01
-1.7029968465264708e-01
-7.0953418018180781e-01
-7.5376948251785270e-01
-7.3385648876540044e-01
-4.9170304973223794e-01
-4.3769486509661026e-01

-2.6641011383243787e+00 9.3367795435260026e-01
2.0146878224597384e-01
-1.6987454338354516e-01
-3.3549431379296490e-01
-3.3419617260933793e-01
-1.9322853570849807e-01
-4.0389286676947966e-02
'''

def test_short_embedding_function_1(potfit):
    potfit.create_param_file()
    potfit.create_potential_file(POTENTIAL.format(1.1, 1.2))
    potfit.create_config_file(repeat_cell=3, seed=42)
    potfit.run()
    assert potfit.has_error()
    assert 'tabulated eqdist' in potfit.stdout
    assert '3 EAM potentials' in potfit.stdout
    assert 'Your embedding function has insufficient sampling points' in potfit.stderr
    assert 'For fixing the gauge degrees of freedom potfit needs to calculate F\'(1.0)!' in potfit.stderr
    assert 'Please include F(1.0) in your potential definition (currently [1.100000,1.200000])' in potfit.stderr

def test_short_embedding_function_2(potfit):
    potfit.create_param_file()
    potfit.create_potential_file(POTENTIAL.format(0.1, 0.2))
    potfit.create_config_file(repeat_cell=3, seed=42)
    potfit.run()
    assert potfit.has_error()
    assert 'tabulated eqdist' in potfit.stdout
    assert '3 EAM potentials' in potfit.stdout
    assert 'Your embedding function has insufficient sampling points' in potfit.stderr
    assert 'For fixing the gauge degrees of freedom potfit needs to calculate F\'(1.0)!' in potfit.stderr
    assert 'Please include F(1.0) in your potential definition (currently [0.100000,0.200000])' in potfit.stderr

def test_short_embedding_function_3(potfit):
    potfit.create_param_file()
    potfit.create_potential_file(POTENTIAL.format(0.1, 1.2))
    potfit.create_config_file(repeat_cell=3, seed=42)
    potfit.run()
    assert potfit.has_no_error()
    assert 'tabulated eqdist' in potfit.stdout
    assert '3 EAM potentials' in potfit.stdout

def test_short_embedding_function_4(potfit):
    potfit.create_param_file()
    potfit.create_potential_file(POTENTIAL.format(1.1, 0.2))
    potfit.create_config_file(repeat_cell=3, seed=42)
    potfit.run()
    assert potfit.has_error()
    assert 'tabulated eqdist' in potfit.stdout
    assert '3 EAM potentials' in potfit.stdout
    assert 'Your embedding function has insufficient sampling points' in potfit.stderr
    assert 'For fixing the gauge degrees of freedom potfit needs to calculate F\'(1.0)!' in potfit.stderr
    assert 'Please include F(1.0) in your potential definition (currently [1.100000,0.200000])' in potfit.stderr
