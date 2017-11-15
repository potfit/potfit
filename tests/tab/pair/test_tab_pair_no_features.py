import pytest

def test_tab_pair_no_features_simple(potfit):
    potfit.create_param_file()
    potfit.create_potential_file('''
#F 3 1
#I 0
#E

1.158337 10.883950976 15

0.00119911169735
7.20484717515e-06
0.000595394899272
0.000936445223073
0.000484087727962
0.000200648587256
0.0000959327293262
0.0000609377489204
0.000049458009818
0.0000457044379028
0.0000444847462971
0.0000440967362355
0.0000439803740159
0.0000439510638516
0.0
''')
    potfit.create_config_file(repeat_cell=3, seed=42)
    potfit.run()
    assert potfit.has_no_error()
    assert 'tabulated eqdist' in potfit.stdout
    assert '1 PAIR potentials' in potfit.stdout
    assert 'Read 1 configuration' in potfit.stdout
    assert 'total of 54 atoms' in potfit.stdout
    assert 'Optimization disabled' in potfit.stdout
    assert 'Potential in format 3 written to file' in potfit.stdout
    assert 'Energy data not written' in potfit.stdout
    assert 'count 163' in potfit.stdout
