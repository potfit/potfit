#!/usr/bin/env python3
################################################################
#
# functions.py:
#   analytic potential functions used in potfit
#
################################################################
#
#   Copyright 2013-2017 -  the potfit development team
#
#   https://www.potfit.net
#
#################################################################
#
#   This file is part of potfit.
#
#   potfit is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   potfit is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with potfit; if not, see <http://www.gnu.org/licenses/>.
#
#################################################################

import random
import sys

if __name__ == '__main__':
    print('Please do not call this script directly!')

class potfit_function(object):
    def __init__(self):
        self.rmin = 0.0
        self.cutoff = 6.0
        self.do_smooth = 0
        self.random = False
        self.params = []

    def set_properties(self, cutoff, do_smooth, random, rmin = 0.0):
        self.cutoff = cutoff
        self.do_smooth = do_smooth
        self.random = random
        self.rmin = rmin
        if self.do_smooth == 1:
            self.params.append(['h', 1, 0.5, 2])

    def set_params(self, params):
        self.params = params

    def write_pot_table(self, output):
        if self.do_smooth != 0:
            output.write('type {}_sc\n'.format(self.__class__.__name__))
        else:
            output.write('type {}\n'.format(self.__class__.__name__))
        output.write('cutoff {}\n'.format(self.cutoff))
        for i in range(len(self.params)):
            output.write('{:12} {}\n'.format(self.params[i][0], self.val_min_max(self.params[i][1], self.params[i][2], self.params[i][3])))
        if self.do_smooth == 2:
            output.write('h!\n')

    def val_min_max(self, val, vmin, vmax):
        if self.random:
            return '{:.2f} {:.2f} {:.2f}'.format(vmin+(vmax-vmin)*random.random(),vmin,vmax)
        else:
            return '{:.2f} {:.2f} {:.2f}'.format(val,vmin,vmax)

    def get_num_params(self):
        if not self.params:
            return -1
        else:
            return len(self.params)

    def get_function_string(self):
        try:
            self.function
        except:
            sys.stderr.write('Error: The analytic function is not defined for {} potentials.\n'.format(self.__class__.__name__))
            sys.exit()
        return self.function

class const(potfit_function):
    def __init__(self):
        super(const, self).__init__()
        self.function = '{0}'
        self.params = [['c', 1, 0, 2]]

class lj(potfit_function):
    def __init__(self):
        super(lj, self).__init__()
        self.function = '4*{0}*(({1}/x)**12-({1}/x)**6)'
        self.params = [
                ['epsilon', 0.1, 0, 1],
                ['sigma', 2.5, 1, 4]]

class morse(potfit_function):
    def __init__(self):
        super(morse, self).__init__()
        self.function = '{0}*((1-exp(-{1}*(x-{2})))**2-1)'
        self.params = [
                ['D_e', 0.1, 0, 1],
                ['a', 2, 1, 5],
                ['r_e', 2.5, 1, 5]]

class power(potfit_function):
    def __init__(self):
        super(power, self).__init__()
        self.function = '{0}*x**{1}'
        self.params = [
                ['alpha', 0.1, 0, 10],
                ['beta', 2, 0, 10]]

class softshell(potfit_function):
    def __init__(self):
        super(softshell, self).__init__()
        self.function = '({0}/x)**{1}'
        self.params = [
                ['alpha', 1, 0.1, 10],
                ['beta', 2, 1, 5]]

class power_decay(potfit_function):
    def __init__(self):
        super(power_decay, self).__init__()
        self.function = '{0}*(1/x)**{1}'
        self.params = [
                ['alpha', 1, 0.1, 10],
                ['beta', 2, 1, 5]]

class exp_decay(potfit_function):
    def __init__(self):
        super(exp_decay, self).__init__()
        self.function = '{0}*exp(-{1}*x)'
        self.params = [
                ['alpha', 4, 1, 20],
                ['beta', 1, 0.5, 5]]

class mexp_decay(potfit_function):
    def __init__(self):
        super(mexp_decay, self).__init__()
        self.function = '{0}*exp(-{1}*(x-{2}))'
        self.params = [
                ['alpha', 0.1, 0, 10],
                ['beta', 0.1, 0, 10],
                ['r_0', 2, 0, 10]]

class exp_plus(potfit_function):
    def __init__(self):
        super(exp_plus, self).__init__()
        self.function = '{0}*exp(-{1}*x)+{2}'
        self.params = [
                ['alpha', 1, 0, 10],
                ['beta', 1, -10, 10],
                ['c', 0, -2, 2]]

class sqrt(potfit_function):
    def __init__(self):
        super(sqrt, self).__init__()
        self.function = '{0}*sqrt(x/{1})'
        self.params = [
                ['alpha', 0.1, 0, 10],
                ['beta', 2, 0, 10]]

class born(potfit_function):
    def __init__(self):
        super(born, self).__init__()
        self.function = '{0}*exp(({4}-x)/{1})-{2}/x**6+{3}/x**8'
        self.params = [
                ['alpha', 0.1, 0, 10],
                ['beta', 2, 0, 10],
                ['gamma', 3, 0, 10],
                ['delta', 2, 0, 10],
                ['r_0', 2.5, 1, 8]]

class harmonic(potfit_function):
    def __init__(self):
        super(harmonic, self).__init__()
        self.function = '{0}*(x-{1})**2'
        self.params = [
                ['alpha', 0.1, 0, 10],
                ['r_0', 2, 0, 8]]

class acosharmonic(potfit_function):
    def __init__(self):
        super(acosharmonic, self).__init__()
        self.function = '{0}*(acos(x)-{1})**2'
        self.params = [
                ['alpha', 0.1, 0, 10],
                ['r_0', 2, 0, 8]]

class eopp(potfit_function):
    def __init__(self):
        super(eopp, self).__init__()
        self.function = '{0}/x**{1}+{2}/x**{3}*cos({4}*x+{5})'
        self.params = [
                ['C_1', 15, 1, 10000],
                ['eta_1', 6, 1, 20],
                ['C_2', 5, -100, 100],
                ['eta_2', 3, 1, 10],
                ['k', 2.5, 0, 6],
                ['phi', 3, 0, 6.3]]

class meopp(potfit_function):
    def __init__(self):
        super(meopp, self).__init__()
        self.function = '{0}/(x-{6})**{1}+{2}/x**{3}*cos({4}*x+{5})'
        self.params = [
                ['C_1', 15, 1, 10000],
                ['eta_1', 6, 1, 20],
                ['C_2', 5, -100, 100],
                ['eta_2', 3, 1, 10],
                ['k', 2.5, 0, 6],
                ['phi', 3, 0, 6.3],
                ['r_0', 0, -3, 3]]

class eopp_exp(potfit_function):
    def __init__(self):
        super(eopp_exp, self).__init__()
        self.function = '{0}*exp(-{1}*x)+{2}/x**{3}*cos({4}*x+{5})'
        self.params = [
                ['C_1', 15, 0.5, 10000],
                ['eta_1', 6, 1, 20],
                ['C_2', 5, -100, 100],
                ['eta_2', 3, 1, 10],
                ['k', 2.5, 0, 6],
                ['phi', 3, 0, 6.3]]

class ms(potfit_function):
    def __init__(self):
        super(ms, self).__init__()
        self.function = '{0}*(exp({1}*(1-x/{2}))-2*exp({1}/2*(1-x/{2})))'
        self.params = [
                ['D_e', 0.1, 0, 1],
                ['a', 2, 1, 5],
                ['r_0', 2.5, 1, 5]]

class strmm(potfit_function):
    def __init__(self):
        super(strmm, self).__init__()
        self.function = '2*{0}*exp(-{1}*(x-{4})/2)-{2}*(1+{3}*(x-{4})*exp(-{3}*(x-{4})))'
        self.params = [
                ['alpha', 1, -10, 10],
                ['beta', 1, -10, 10],
                ['gamma', 1, -10, 10],
                ['delta', 1, -10, 10],
                ['r_0', 1, -10, 10]]

class double_morse(potfit_function):
    def __init__(self):
        super(double_morse, self).__init__()
        self.function = '{0}*(exp(-2*{1}*(x-{2}))-2*exp(-{1}*(x-{2})))+{3}*(exp(-2*{4}*(x-{5}))-2*exp(-{4}*(x-{5})))+{6}'
        self.params = [
                ['E_1', 1, -10, 10],
                ['alpha_1', 1, -10, 10],
                ['r_0', 1, -10, 10],
                ['E_2', 1, -10, 10],
                ['alpha_2', 1, -10, 10],
                ['r_1', 1, -10, 10]
                ['delta', 1, -10, 10]]

class double_exp(potfit_function):
    def __init__(self):
        super(double_exp, self).__init__()
        self.function = '{0}*exp(-{1}*(x-{2})**2)+exp(-{3}*(x-{4}))'
        self.params = [
                ['a', 1, -10, 10],
                ['beta_1', 1, -10, 10],
                ['r_0', 1, -10, 10],
                ['beta_2', 1, -10, 10],
                ['r_1', 1, -10, 10]]

class poly_5(potfit_function):
    def __init__(self):
        super(poly_5, self).__init__()
        self.function = '{0}+0.5*{1}*(x-1)**2+{2}*(x-1)**3+{3}*(r-1)**4+{5}*(r-1)**5'
        self.params = [
                ['F_0', 1, -10, 10],
                ['F_2', 1, -10, 10],
                ['q_1', 1, -10, 10],
                ['q_2', 1, -10, 10],
                ['q_3', 1, -10, 10]]

class buck(potfit_function):
    def __init__(self):
        super(buck, self).__init__()
        self.function = '{0}*exp(-x/{1})-{2}*({1}/x)**6'
        self.params = [
                ['alpha', 1, -10, 10],
                ['beta', 1, -10, 10],
                ['gamma', 1, -10, 10]]

class kawamura(potfit_function):
    def __init__(self):
        super(kawamura, self).__init__()
        self.function = '{0}*{1}/x+{2}*({5}+{6})*exp(({3}+{4}-x)/({5}+{6}))-{7}*{8}/x**6'
        self.params = [
                ['z_1', 1, -10, 10],
                ['z_2', 1, -10, 10],
                ['f_0', 1, -10, 10],
                ['a_1', 1, -10, 10],
                ['a_2', 1, -10, 10],
                ['b_1', 1, -10, 10],
                ['b_2', 1, -10, 10],
                ['c_1', 1, -10, 10],
                ['c_2', 1, -10, 10]]

class kawamura_mix(potfit_function):
    def __init__(self):
        super(kawamura_mix, self).__init__()
        self.function = '{0}*{1}/x+{2}*({5}+{6})*exp(({3}+{4}-x)/({5}+{6}))-{7}*{8}/x**6+{2}*{9}*(exp(-2*{10}*(x-{11}))-2*exp(-{10}*(x-{11})))'
        self.params = [
                ['z_1', 1, -10, 10],
                ['z_2', 1, -10, 10],
                ['f_0', 1, -10, 10],
                ['a_1', 1, -10, 10],
                ['a_2', 1, -10, 10],
                ['b_1', 1, -10, 10],
                ['b_2', 1, -10, 10],
                ['c_1', 1, -10, 10],
                ['c_2', 1, -10, 10],
                ['D', 1, -10, 10],
                ['beta', 1, -10, 10],
                ['r_0', 1, -10, 10]]

class mishin(potfit_function):
    def __init__(self):
        super(mishin, self).__init__()
        self.function = '{0}*(x-{3})**{4}*exp(-{5}*(x-{3}))*(1+{1}*exp(-{5}*(x-{3})))+{2}'
        self.params = [
                ['A_0', 1, -10, 10],
                ['B_0', 1, -10, 10],
                ['C_0', 1, -10, 10],
                ['r_0', 1, -10, 10],
                ['y', 1, -10, 10],
                ['gamma', 1, -10, 10]]

class gen_lj(potfit_function):
    def __init__(self):
        super(gen_lj, self).__init__()
        self.function = '{0}/({2}-{1})*({2}/(r/{3})**{1}-{1}/(x/{3})**{2})+{4}'
        self.params = [
                ['V_0', 1, -10, 10],
                ['b_1', 1, -10, 10],
                ['b_2', 1, -10, 10],
                ['r_1', 1, -10, 10],
                ['delta', 1, -10, 10]]

class gljm(potfit_function):
    def __init__(self):
        super(gljm, self).__init__()
        self.function = '{0}/({2}-{1})*({2}/(r/{3})**{1}-{1}/(x/{3})**{2})+{4}+{5}*({6}*(x-{9})**{10}*exp(-{11}*(x-{9}))*(1+{7}*exp(-{11}*(x-{9})))+{8})'
        self.params = [
                ['V_0', 1, -10, 10],
                ['b_1', 1, -10, 10],
                ['b_2', 1, -10, 10],
                ['r_1', 1, -10, 10],
                ['delta', 1, -10, 10],
                ['A_0', 1, -10, 10],
                ['B_0', 1, -10, 10],
                ['C_0', 1, -10, 10],
                ['r_0', 1, -10, 10],
                ['y', 1, -10, 10],
                ['gamma', 1, -10, 10]]

class vas(potfit_function):
    def __init__(self):
        super(vas, self).__init__()
        self.function = 'exp({0}/(x-{1}))'
        self.params = [
                ['alpha', 1, -10, 10],
                ['beta', 1, -10, 10]]

class vpair(potfit_function):
    def __init__(self):
        super(vpair, self).__init__()
        self.function = '14.4*({0}/x**{1}-({4}*{3}*{3}+{5}*{2}*{2})/x**4*exp(-x/{6}))'
        self.params = [
                ['alpha', 1, -10, 10],
                ['beta', 1, -10, 10],
                ['gamma', 1, -10, 10],
                ['delta', 1, -10, 10],
                ['a', 1, -10, 10],
                ['b', 1, -10, 10],
                ['c', 1, -10, 10]]

class universal(potfit_function):
    def __init__(self):
        super(universal, self).__init__()
        self.function = '{0}*({2}/({2}-{1})*x**{1}-{1}/({2}-{1})*x**{2})+{3}*x'
        self.params = [
                ['F_0', 1, -10, 10],
                ['m', 1, 0, 20],
                ['n', 2, 0, 20],
                ['F_1', 0, -1, 1]]

class bjs(potfit_function):
    def __init__(self):
        super(bjs, self).__init__()
        self.function = '{0}*(1-{1}*log(x))*x**{1}+{2}'
        self.params = [
                ['F_0', -1, -10, 0],
                ['gamma', 2, 0.1, 2],
                ['F_1', 2, 1, 5]]

class parabola(potfit_function):
    def __init__(self):
        super(parabola, self).__init__()
        self.function = '{0}*x**2+{1}*x+{2}'
        self.params = [
                ['alpha', 1, -10, 10],
                ['beta', 1, -10, 10],
                ['gamma', 1, -10, 10]]

class csw(potfit_function):
    def __init__(self):
        super(csw, self).__init__()
        self.function = '(1+{0}*cos({2}*x)+{1}*sin({2}*x))/x**{3}'
        self.params = [
                ['a_1', 0.2, -2, 2],
                ['a_2', 0.2, -2, 2],
                ['alpha', 2, 1, 6],
                ['beta', 3, 0.5, 5]]

class csw2(potfit_function):
    def __init__(self):
        super(csw2, self).__init__()
        self.function = '(1+{0}*cos({1}*x+{2}))/x**{3}'
        self.params = [
                ['a', 0.2, -2, 2],
                ['alpha', 2, 1, 6],
                ['phi', 0, 0, 6.3],
                ['beta', 3, 0.5, 5]]

class tersoff_pot(potfit_function):
    def __init__(self):
        super(tersoff_pot, self).__init__()
        self.params = [
                ['A', 1, -10, 10],
                ['B', 1, -10, 10],
                ['lambda', 1, -10, 10],
                ['mu', 1, -10, 10],
                ['gamma', 1, -10, 10],
                ['n', 1, -10, 10],
                ['c', 1, -10, 10],
                ['d', 1, -10, 10],
                ['h', 1, -10, 10],
                ['S', 1, -10, 10],
                ['R', 1, -10, 10]]

class tersoff_mix(potfit_function):
    def __init__(self):
        super(tersoff_mix, self).__init__()
        self.params = [
                ['chi', 1, -10, 10],
                ['omega', 1, -10, 10]]

class stiweb_2(potfit_function):
    def __init__(self, cutoff):
        super(stiweb_2, self).__init__()
        self.params = [
                ['A', 1, -10, 10],
                ['B', 1, -10, 10],
                ['p', 1, -10, 10],
                ['q', 1, -10, 10],
                ['delta', 1, -10, 10],
                ['a', 1, 0, cutoff]]

class stiweb_3(potfit_function):
    def __init__(self, cutoff):
        super(stiweb_3, self).__init__()
        self.params = [
                ['gamma', 1, -10, 10],
                ['b', 1, 0, cutoff]]
