function out = model
%
% test.m
%
% Model exported on Feb 21 2020, 14:55 by COMSOL 5.4.0.246.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('D:\Christian Haffner\Simulations\NLOSPG\Template');

model.label('OEO_Slot.mph');

model.param.set('wSim', '6 [um]');
model.param.set('hSubstrate', '4 [um]');
model.param.set('hWG_top', '2.2e-07 [m]');
model.param.set('hWG_bot', '160 [nm]');
model.param.set('wWG', '500 [nm]');
model.param.set('hOrganic', '2e-07 [m]');
model.param.set('hBuffer', '0 [m]');
model.param.set('hOEO', '50 [nm]');
model.param.set('hcladding', '2000 [nm]');
model.param.set('r33', '137 [pm/V]');
model.param.set('lmin', '1100 [nm]');
model.param.set('lmax', '2300 [nm]');
model.param.set('Nl', '25  ');
model.param.set('f_rf', '9000000000 [Hz]');
model.param.set('n_start', '3.475  ');
model.param.set('nr_modes', '3  ');
model.param.set('V_bias', '1 [V]');
model.param.set('ASiO2', '0.69617  ');
model.param.set('BSiO2', '0.40794  ');
model.param.set('CSiO2', '0.89748  ');
model.param.set('l0SiO2', '0.068404 [um]');
model.param.set('l1SiO2', '0.11624 [um]');
model.param.set('l2SiO2', '9.8961 [um]');
model.param.set('lp', '0.14586 [um]');
model.param.set('Gamma', '71428571428571.42 [Hz]');
model.param.set('ASiNx', '3.0249  ');
model.param.set('BSiNx', '40314  ');
model.param.set('l0SiNx', '0.13534 [um]');
model.param.set('l1SiNx', '1239.842 [um]');
model.param.set('AAl2O3', '1.2038 [um]');
model.param.set('BAl2O3', '1.0583 [um]');
model.param.set('CAl2O3', '5.2807 [um]');
model.param.set('l0Al2O3', '0.061448 [um]');
model.param.set('l1Al2O3', '0.1107 [um]');
model.param.set('l2Al2O3', '17.9266 [um]');
model.param.set('e', '1.602e-19 [C]');
model.param.set('hbar', '1.0546e-34 [J*s]');
model.param.set('eps0', '8.854e-12 [F/m]');
model.param.set('mu0', '1.2566e-06 [H/m]');
model.param.set('c_0', '300000000 [m/s]');
model.param.set('material_passive_photonic', '300000000 [m/s]');
model.param.set('wl', '1.55e-06 [m]');
model.param.set('plasmonic_mesh', '150  ');
model.param.set('omega', '2*pi*c_const/wl ');
model.param.set('gamma', '0.001*omega ');
model.param.set('sigmaxx', '-1i*e^2/(4*pi*hbar)*log((-omega + 1i*2*gamma)/(+omega + 1i*2*gamma)) ');
model.param.set('sigmayy', '0 [S]');
model.param.set('sigmazz', '-1i*e^2/(4*pi*hbar)*log((-omega + 1i*2*gamma)/(+omega + 1i*2*gamma)) ');

model.component.create('comp1', true);

model.func.create('Graphene', 'Analytic');
model.func('Graphene').active(true);
model.func('Graphene').set('expr', 'imag(sqrt(1-(wl*1e-9)^2/lp^2+1i*Gamma*(wl*1e-9)^3/(2*pi*c_const*lp^2)))');
model.func('Graphene').set('args', 'wl');
model.func('Graphene').set('argunit', 'wl');
model.func('Graphene').set('plotargs', {'wl' '500' '2000'});
model.func.create('Substrate', 'Interpolation');
model.func('Substrate').set('funcs', {'epsSubstrate_re' '1'; 'epsSubstrate_im' '2'});
model.func('Substrate').set('source', 'file');
model.func('Substrate').set('filename', 'D:\Christian Haffner\Simulations\NLOSPG\DataIn\eps_SiO2_Lemarchand_2013_250nm__2500nm.txt');
model.func('Substrate').set('extrap', 'linear');
model.func.create('PhotonicWG', 'Interpolation');
model.func('PhotonicWG').set('funcs', {'epsPhotonicWG_re' '1'; 'epsPhotonicWG_im' '2'});
model.func('PhotonicWG').set('source', 'file');
model.func('PhotonicWG').set('filename', 'D:\Christian Haffner\Simulations\NLOSPG\DataIn\eps_Si_Green_2008.txt');
model.func('PhotonicWG').set('extrap', 'linear');
model.func.create('OEO', 'Interpolation');
model.func('OEO').set('funcs', {'epsOEO_re' '1'; 'epsOEO_im' '2'});
model.func('OEO').set('source', 'file');
model.func('OEO').set('filename', 'D:\Christian Haffner\Simulations\NLOSPG\DataIn\eps_HD_BB_OH_UniWashington.txt');
model.func('OEO').set('extrap', 'linear');
model.func.create('Metal_1', 'Interpolation');
model.func('Metal_1').set('funcs', {'epsMetal_1_re' '1'; 'epsMetal_1_im' '2'});
model.func('Metal_1').set('source', 'file');
model.func('Metal_1').set('filename', 'D:\Christian Haffner\Simulations\NLOSPG\DataIn\eps_Au_IEF_1604.txt');
model.func('Metal_1').set('extrap', 'linear');
model.func.create('Metal_2', 'Interpolation');
model.func('Metal_2').set('funcs', {'epsMetal_2_re' '1'; 'epsMetal_2_im' '2'});
model.func('Metal_2').set('source', 'file');
model.func('Metal_2').set('filename', 'D:\Christian Haffner\Simulations\NLOSPG\DataIn\eps_Cu_McPeak.txt');
model.func('Metal_2').set('extrap', 'linear');

model.component('comp1').geom.create('geom1', 2);
model.component('comp1').geom('geom1').selection.create('Substrate', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('Substrate').label('Substrate');
model.component('comp1').geom('geom1').selection.create('PhotonicWG', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('PhotonicWG').label('PhotonicWG');
model.component('comp1').geom('geom1').selection.create('OEO', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('OEO').label('OEO');
model.component('comp1').geom('geom1').selection.create('Metal_1', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('Metal_1').label('Metal_1');
model.component('comp1').geom('geom1').selection.create('Metal_2', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('Metal_2').label('Metal_2');
model.component('comp1').geom('geom1').selection.create('Graphene', 'CumulativeSelection');
model.component('comp1').geom('geom1').selection('Graphene').label('Graphene');
model.component('comp1').geom('geom1').create('r_ThermalOxide', 'Rectangle');
model.component('comp1').geom('geom1').feature('r_ThermalOxide').label('Thermal Oxide');
model.component('comp1').geom('geom1').feature('r_ThermalOxide').set('base', 'center');
model.component('comp1').geom('geom1').feature('r_ThermalOxide').set('pos', {'0' '-hSubstrate/2'});
model.component('comp1').geom('geom1').feature('r_ThermalOxide').set('size', {'wSim' 'hSubstrate'});
model.component('comp1').geom('geom1').feature('r_ThermalOxide').set('contributeto', 'Substrate');
model.component('comp1').geom('geom1').create('r_cladding', 'Rectangle');
model.component('comp1').geom('geom1').feature('r_cladding').label('cladding');
model.component('comp1').geom('geom1').feature('r_cladding').set('base', 'center');
model.component('comp1').geom('geom1').feature('r_cladding').set('pos', {'0' 'hcladding/2+hOEO'});
model.component('comp1').geom('geom1').feature('r_cladding').set('size', {'wSim' 'hcladding'});
model.component('comp1').geom('geom1').feature('r_cladding').set('contributeto', 'Substrate');
model.component('comp1').geom('geom1').create('r_photonic_wg_bot', 'Rectangle');
model.component('comp1').geom('geom1').feature('r_photonic_wg_bot').label('photonic_wg_bot');
model.component('comp1').geom('geom1').feature('r_photonic_wg_bot').set('base', 'center');
model.component('comp1').geom('geom1').feature('r_photonic_wg_bot').set('pos', {'0' '-hWG_bot/2'});
model.component('comp1').geom('geom1').feature('r_photonic_wg_bot').set('size', {'wWG' 'hWG_bot'});
model.component('comp1').geom('geom1').feature('r_photonic_wg_bot').set('contributeto', 'PhotonicWG');
model.component('comp1').geom('geom1').create('r_photonic_wg_top', 'Rectangle');
model.component('comp1').geom('geom1').feature('r_photonic_wg_top').label('photonic_wg_top');
model.component('comp1').geom('geom1').feature('r_photonic_wg_top').set('base', 'center');
model.component('comp1').geom('geom1').feature('r_photonic_wg_top').set('pos', {'0' 'hWG_top/2 + hOEO'});
model.component('comp1').geom('geom1').feature('r_photonic_wg_top').set('size', {'wWG' 'hWG_top'});
model.component('comp1').geom('geom1').feature('r_photonic_wg_top').set('contributeto', 'PhotonicWG');
model.component('comp1').geom('geom1').create('r_OEO_slot', 'Rectangle');
model.component('comp1').geom('geom1').feature('r_OEO_slot').label('OEO_slot');
model.component('comp1').geom('geom1').feature('r_OEO_slot').set('base', 'center');
model.component('comp1').geom('geom1').feature('r_OEO_slot').set('pos', {'0' 'hOEO/2'});
model.component('comp1').geom('geom1').feature('r_OEO_slot').set('size', {'wSim' 'hOEO'});
model.component('comp1').geom('geom1').feature('r_OEO_slot').set('contributeto', 'OEO');
model.component('comp1').geom('geom1').create('poly_graphene_bot', 'Polygon');
model.component('comp1').geom('geom1').feature('poly_graphene_bot').set('source', 'vectors');
model.component('comp1').geom('geom1').feature('poly_graphene_bot').set('type', 'solid');
model.component('comp1').geom('geom1').feature('poly_graphene_bot').set('x', '-wSim/2 wWG/2');
model.component('comp1').geom('geom1').feature('poly_graphene_bot').set('y', '0 0');
model.component('comp1').geom('geom1').feature('poly_graphene_bot').set('contributeto', 'Graphene');
model.component('comp1').geom('geom1').create('poly_graphene_top', 'Polygon');
model.component('comp1').geom('geom1').feature('poly_graphene_top').set('source', 'vectors');
model.component('comp1').geom('geom1').feature('poly_graphene_top').set('type', 'solid');
model.component('comp1').geom('geom1').feature('poly_graphene_top').set('x', '-wWG/2 wSim/2');
model.component('comp1').geom('geom1').feature('poly_graphene_top').set('y', 'hOEO hOEO');
model.component('comp1').geom('geom1').feature('poly_graphene_top').set('contributeto', 'Graphene');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run;

model.component('comp1').material.create('matSubstrate', 'Common');
model.component('comp1').material('matSubstrate').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('matSubstrate').label('Substrate');
model.component('comp1').material('matSubstrate').propertyGroup('RefractiveIndex').set('n', []);
model.component('comp1').material('matSubstrate').propertyGroup('RefractiveIndex').set('ki', []);
model.component('comp1').material('matSubstrate').propertyGroup('RefractiveIndex').set('n', {'real(sqrt(epsSubstrate_re(wl*1e9)+i*epsSubstrate_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsSubstrate_re(wl*1e9)+i*epsSubstrate_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsSubstrate_re(wl*1e9)+i*epsSubstrate_im(wl*1e9)))'});
model.component('comp1').material('matSubstrate').propertyGroup('RefractiveIndex').set('ki', {'imag(sqrt(epsSubstrate_re(wl*1e9)+i*epsSubstrate_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsSubstrate_re(wl*1e9)+i*epsSubstrate_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsSubstrate_re(wl*1e9)+i*epsSubstrate_im(wl*1e9)))'});
model.component('comp1').material('matSubstrate').selection.set([1 3 4 5]);
model.component('comp1').material.create('matPhotonicWG', 'Common');
model.component('comp1').material('matPhotonicWG').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('matPhotonicWG').label('PhotonicWG');
model.component('comp1').material('matPhotonicWG').propertyGroup('RefractiveIndex').set('n', []);
model.component('comp1').material('matPhotonicWG').propertyGroup('RefractiveIndex').set('ki', []);
model.component('comp1').material('matPhotonicWG').propertyGroup('RefractiveIndex').set('n', {'real(sqrt(epsPhotonicWG_re(wl*1e9)+i*epsPhotonicWG_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsPhotonicWG_re(wl*1e9)+i*epsPhotonicWG_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsPhotonicWG_re(wl*1e9)+i*epsPhotonicWG_im(wl*1e9)))'});
model.component('comp1').material('matPhotonicWG').propertyGroup('RefractiveIndex').set('ki', {'imag(sqrt(epsPhotonicWG_re(wl*1e9)+i*epsPhotonicWG_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsPhotonicWG_re(wl*1e9)+i*epsPhotonicWG_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsPhotonicWG_re(wl*1e9)+i*epsPhotonicWG_im(wl*1e9)))'});
model.component('comp1').material('matPhotonicWG').selection.set([4 5]);
model.component('comp1').material.create('matOEO', 'Common');
model.component('comp1').material('matOEO').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('matOEO').label('OEO');
model.component('comp1').material('matOEO').propertyGroup('RefractiveIndex').set('n', []);
model.component('comp1').material('matOEO').propertyGroup('RefractiveIndex').set('ki', []);
model.component('comp1').material('matOEO').propertyGroup('RefractiveIndex').set('n', {'real(sqrt(epsOEO_re(wl*1e9)+i*epsOEO_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsOEO_re(wl*1e9)+i*epsOEO_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsOEO_re(wl*1e9)+i*epsOEO_im(wl*1e9)))'});
model.component('comp1').material('matOEO').propertyGroup('RefractiveIndex').set('ki', {'imag(sqrt(epsOEO_re(wl*1e9)+i*epsOEO_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsOEO_re(wl*1e9)+i*epsOEO_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsOEO_re(wl*1e9)+i*epsOEO_im(wl*1e9)))'});
model.component('comp1').material('matOEO').selection.set([2]);
model.component('comp1').material('matOEO').propertyGroup('def').set('relpermittivity', '5.6');
model.component('comp1').material.create('matMetal_1', 'Common');
model.component('comp1').material('matMetal_1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('matMetal_1').label('Metal_1');
model.component('comp1').material('matMetal_1').propertyGroup('RefractiveIndex').set('n', []);
model.component('comp1').material('matMetal_1').propertyGroup('RefractiveIndex').set('ki', []);
model.component('comp1').material('matMetal_1').propertyGroup('RefractiveIndex').set('n', {'real(sqrt(epsMetal_1_re(wl*1e9)+i*epsMetal_1_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsMetal_1_re(wl*1e9)+i*epsMetal_1_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsMetal_1_re(wl*1e9)+i*epsMetal_1_im(wl*1e9)))'});
model.component('comp1').material('matMetal_1').propertyGroup('RefractiveIndex').set('ki', {'imag(sqrt(epsMetal_1_re(wl*1e9)+i*epsMetal_1_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsMetal_1_re(wl*1e9)+i*epsMetal_1_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsMetal_1_re(wl*1e9)+i*epsMetal_1_im(wl*1e9)))'});
model.component('comp1').material('matMetal_1').selection.set([]);
model.component('comp1').material.create('matMetal_2', 'Common');
model.component('comp1').material('matMetal_2').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('matMetal_2').label('Metal_2');
model.component('comp1').material('matMetal_2').propertyGroup('RefractiveIndex').set('n', []);
model.component('comp1').material('matMetal_2').propertyGroup('RefractiveIndex').set('ki', []);
model.component('comp1').material('matMetal_2').propertyGroup('RefractiveIndex').set('n', {'real(sqrt(epsMetal_2_re(wl*1e9)+i*epsMetal_2_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsMetal_2_re(wl*1e9)+i*epsMetal_2_im(wl*1e9)))' '0' '0' '0' 'real(sqrt(epsMetal_2_re(wl*1e9)+i*epsMetal_2_im(wl*1e9)))'});
model.component('comp1').material('matMetal_2').propertyGroup('RefractiveIndex').set('ki', {'imag(sqrt(epsMetal_2_re(wl*1e9)+i*epsMetal_2_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsMetal_2_re(wl*1e9)+i*epsMetal_2_im(wl*1e9)))' '0' '0' '0' 'imag(sqrt(epsMetal_2_re(wl*1e9)+i*epsMetal_2_im(wl*1e9)))'});
model.component('comp1').material('matMetal_2').selection.set([]);
model.component('comp1').material.create('matGraphene', 'Common');
model.component('comp1').material('matGraphene').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('matGraphene').label('Graphene');
model.component('comp1').material('matGraphene').propertyGroup('RefractiveIndex').set('n', []);
model.component('comp1').material('matGraphene').propertyGroup('RefractiveIndex').set('ki', []);
model.component('comp1').material('matGraphene').propertyGroup('RefractiveIndex').set('n', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.component('comp1').material('matGraphene').propertyGroup('RefractiveIndex').set('ki', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('matGraphene').selection.set([]);

model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').create('ftriMetal_2', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').label('ftriMetal_2');
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').set('xscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').set('yscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').feature('size1').set('hmax', 'wl/plasmonic_mesh');
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftriMetal_2').selection.set([]);
model.component('comp1').mesh('mesh1').create('ftriMetal_1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').label('ftriMetal_1');
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').set('xscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').set('yscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').feature('size1').set('hmax', 'wl/plasmonic_mesh');
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftriMetal_1').selection.set([]);
model.component('comp1').mesh('mesh1').create('ftriOEO', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftriOEO').label('ftriOEO');
model.component('comp1').mesh('mesh1').feature('ftriOEO').set('xscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriOEO').set('yscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriOEO').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftriOEO').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftriOEO').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftriOEO').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftriOEO').feature('size1').set('hmax', 'wl/1.7701/8');
model.component('comp1').mesh('mesh1').feature('ftriOEO').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftriOEO').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftriOEO').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftriOEO').selection.set([2]);
model.component('comp1').mesh('mesh1').create('ftriPhotonicWG', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').label('ftriPhotonicWG');
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').set('xscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').set('yscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').feature('size1').set('hmax', 'wl/3.475/8');
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftriPhotonicWG').selection.set([4 5]);
model.component('comp1').mesh('mesh1').create('ftriSubstrate', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').label('ftriSubstrate');
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').set('xscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').set('yscale', 1);
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').feature('size1').set('hmax', 'wl/1.4636/8');
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftriSubstrate').selection.set([1 3]);
model.component('comp1').mesh('mesh1').run;

model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
model.component('comp1').physics('ewfd').create('sctr1', 'Scattering', 1);

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('ewfd').feature('sctr1').selection.set([1 2 4 5 6 7 10 12 15 17 18 20 3 4 6 10 12 15 17 19]);
model.component('comp1').physics('ewfd').create('scu1', 'SurfaceCurrent', 1);
model.component('comp1').physics('ewfd').feature('scu1').selection.set([4 10 12 17]);
model.component('comp1').physics('ewfd').feature('scu1').set('Js0', {'sigmaxx*ewfd.Ex' 'sigmayy*ewfd.Ey' 'sigmazz*ewfd.Ez'});
model.component('comp1').physics.create('es', 'Electrostatics', 'geom1');
model.component('comp1').physics('es').selection.set([2]);
model.component('comp1').physics('es').create('gnd1', 'Ground', 1);
model.component('comp1').physics('es').feature('gnd1').selection.set([12 17]);
model.component('comp1').physics('es').create('pot1', 'ElectricPotential', 1);
model.component('comp1').physics('es').feature('pot1').selection.set([4 10]);
model.component('comp1').physics('es').feature('pot1').set('V0', 'V_bias');

model.study.create('std1');
model.study('std1').setGenConv(true);
model.study('std1').create('stat', 'Stationary');
model.study('std1').feature('stat').activate('ewfd', false);
model.study('std1').feature('stat').activate('es', true);
model.study('std1').feature('stat').setIndex('activate', false, 1);
model.study('std1').create('mode', 'ModeAnalysis');
model.study('std1').feature('mode').set('ngen', '5');
model.study('std1').feature('mode').activate('ewfd', true);
model.study('std1').feature('mode').set('modeFreq', 'c_const/wl');
model.study('std1').feature('mode').set('shiftactive', true);
model.study('std1').feature('mode').set('shift', 'n_start');
model.study('std1').feature('mode').set('neigsactive', true);
model.study('std1').feature('mode').set('neigs', 3);

model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');
model.study('std1').feature('mode').set('notlistsolnum', 1);
model.study('std1').feature('mode').set('notsolnum', 'auto');
model.study('std1').feature('mode').set('listsolnum', 1);
model.study('std1').feature('mode').set('solnum', 'auto');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').create('su1', 'StoreSolution');
model.sol('sol1').create('st2', 'StudyStep');
model.sol('sol1').feature('st2').set('study', 'std1');
model.sol('sol1').feature('st2').set('studystep', 'mode');
model.sol('sol1').create('v2', 'Variables');
model.sol('sol1').feature('v2').set('initmethod', 'sol');
model.sol('sol1').feature('v2').set('initsol', 'sol1');
model.sol('sol1').feature('v2').set('initsoluse', 'su1');
model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
model.sol('sol1').feature('v2').set('notsol', 'sol1');
model.sol('sol1').feature('v2').set('control', 'mode');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').set('neigs', 6);
model.sol('sol1').feature('e1').set('shift', '1');
model.sol('sol1').feature('e1').set('control', 'mode');
model.sol('sol1').feature('e1').set('linpmethod', 'sol');
model.sol('sol1').feature('e1').set('linpsol', 'sol1');
model.sol('sol1').feature('e1').set('linpsoluse', 'su1');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol('sol1').feature('v2').set('notsolnum', 'auto');
model.sol('sol1').feature('v2').set('notsolvertype', 'solnum');
model.sol('sol1').feature('v2').set('notlistsolnum', '1');
model.sol('sol1').feature('v2').set('notsolnum', 'auto');
model.sol('sol1').feature('v2').set('notlistsolnum', '1');
model.sol('sol1').feature('v2').set('notsolnum', 'auto');
model.sol('sol1').feature('v2').set('control', 'mode');
model.sol('sol1').attach('std1');

model.label('testsets22.mph');

model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').label('Electric Field (ewfd)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').set('data', 'dset1');
model.result('pg1').feature.create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('data', 'parent');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg2').label('Electric Potential (es)');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').set('data', 'dset1');
model.result('pg2').feature.create('surf1', 'Surface');
model.result('pg2').feature('surf1').set('expr', 'V');
model.result('pg2').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg2').feature('surf1').set('data', 'parent');

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg1').set('looplevel', [2]);
model.result('pg1').run;
model.result('pg1').set('looplevel', [3]);
model.result('pg1').run;
model.result('pg1').set('looplevel', [1]);
model.result('pg1').run;

model.study('std1').create('param', 'Parametric');
model.study('std1').feature.remove('param');

model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').set('data', 'dset2');
model.result('pg2').run;
model.result('pg1').run;
model.result('pg1').run;
model.result.numerical.create('int1', 'IntSurface');
model.result.numerical('int1').set('intvolume', true);
model.result.numerical('int1').setIndex('expr', 'sqrt(eps0*(ewfd.nxx^2*(abs(ewfd.Ex)^2 + ewfd.nyy^2*abs(ewfd.Ey)^2 + ewfd.nzz^2*abs(ewfd.Ez)^2)))', 0);
model.result.numerical('int1').selection.set([1 2 3 4 5]);
model.result.numerical('int1').setIndex('expr', '(eps0*(ewfd.nxx^2*(abs(ewfd.Ex)^2 + ewfd.nyy^2*abs(ewfd.Ey)^2 + ewfd.nzz^2*abs(ewfd.Ez)^2)))', 0);
model.result.numerical('int1').setIndex('expr', 'sqrt(eps0*(ewfd.nxx^2*(abs(ewfd.Ex)^2 + ewfd.nyy^2*abs(ewfd.Ey)^2 + ewfd.nzz^2*abs(ewfd.Ez)^2)))', 0);
model.result.numerical('int1').setIndex('expr', '(eps0*(ewfd.nxx^2*(abs(ewfd.Ex)^2 + ewfd.nyy^2*abs(ewfd.Ey)^2 + ewfd.nzz^2*abs(ewfd.Ez)^2)))', 0);
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('gev1').set('expr', {'es.zref'});
model.result.numerical('gev1').set('descr', {'Reference impedance'});
model.result.numerical('gev1').set('unit', {['ohm' ]});
model.result.table.create('tbl1', 'Table');
model.result.table('tbl1').comments('Global Evaluation 1 (es.zref)');
model.result.numerical('gev1').set('table', 'tbl1');
model.result.numerical('gev1').setResult;
model.result.numerical.remove('gev1');
model.result.numerical.create('gmev1', 'EvalGlobalMatrix');
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('gev1').setIndex('expr', '(hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO)/V_bias^2)', 0);
model.result.table.create('tbl2', 'Table');
model.result.table('tbl2').comments('Global Evaluation 1 ((hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO)/V_bias^2))');
model.result.numerical('gev1').set('table', 'tbl2');
model.result.numerical('gev1').set('data', 'dset2');
model.result.table.create('tbl3', 'Table');
model.result.table('tbl3').comments('Global Evaluation 1 ((hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO)/V_bias^2))');
model.result.numerical('gev1').set('table', 'tbl3');
model.result.table.create('tbl4', 'Table');
model.result.table('tbl4').comments('Global Evaluation 1 ((hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO)/V_bias^2))');
model.result.numerical('gev1').set('table', 'tbl4');
model.result.numerical('int1').setIndex('expr', '(hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO)/V_bias^2)', 1);
model.result.numerical('int1').setIndex('expr', '(hbar*omega*2/(eps0*es.epsilonryy*wWG^2/hOEO)/V_bias^2)', 1);
model.result.numerical('int1').setIndex('expr', 'hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO)/V_bias^2', 1);
model.result.numerical('int1').setIndex('expr', 'hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO*V_bias^2)', 1);
model.result.numerical('int1').setIndex('expr', 'omega', 2);
model.result.numerical('int1').setIndex('expr', '', 2);
model.result.numerical('int1').setIndex('expr', '2*pi*c_const/wl', 2);
model.result.numerical('int1').setIndex('expr', 'c_const/wl', 2);
model.result.numerical('int1').setIndex('expr', 'c_const', 2);

model.param.set('omega', '2*pi*c_0/wl');

model.result.numerical('int1').setIndex('expr', '2*pi*c_0/wl', 2);
model.result.numerical('gev1').setIndex('expr', '2*pi*c_0/wl', 1);
model.result.numerical('gev1').setIndex('expr', '', 1);
model.result.numerical('gev1').setIndex('expr', 'hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO*V_bias^2)', 1);
model.result.numerical('gev1').setIndex('expr', 'omega', 2);
model.result.numerical('gev1').setIndex('expr', 'hbar*omega', 2);
model.result.numerical('gev1').setIndex('expr', 'hbar*omega/eps0', 2);
model.result.numerical('gev1').setIndex('expr', 'hbar*omega/eps0/V_bias', 2);
model.result.numerical('gev1').setIndex('expr', 'hbar*omega/eps0/V_bias^2', 2);
model.result.numerical('gev1').setIndex('expr', 'eps0/V_bias^2', 2);
model.result.numerical('gev1').setIndex('expr', '1/eps0*V_bias^2', 2);
model.result.numerical('gev1').setIndex('expr', 'eps0*V_bias^2', 2);
model.result.numerical('gev1').setIndex('expr', 'eps0*es.epsilonryy*wWG/hOEO*V_bias^2)/(hbar*omega*2)', 1);
model.result.numerical('gev1').setIndex('expr', '(eps0*es.epsilonryy*wWG/hOEO*V_bias^2)/(hbar*omega*2)', 1);
model.result.numerical('gev1').setIndex('expr', '', 2);
model.result.numerical('gev1').setIndex('expr', '(eps0*es.epsilonryy*wWG/hOEO*V_bias^2)', 1);
model.result.table.create('tbl5', 'Table');
model.result.table('tbl5').comments('Global Evaluation 1 ((hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO)/V_bias^2), (eps0*es.epsilonryy*wWG/hOEO*V_bias^2), )');
model.result.numerical('gev1').set('table', 'tbl5');
model.result.table.create('tbl6', 'Table');
model.result.table('tbl6').comments('Global Evaluation 1 ((hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO)/V_bias^2), (eps0*es.epsilonryy*wWG/hOEO*V_bias^2), )');
model.result.numerical('gev1').set('table', 'tbl6');
model.result.numerical('int1').setIndex('expr', '(eps0*es.epsilonryy*wWG/hOEO*V_bias^2)', 2);
model.result.table.create('tbl7', 'Table');
model.result.table('tbl7').comments('Surface Integration 1 ((eps0*(ewfd.nxx^2*(abs(ewfd.Ex)^2 + ewfd.nyy^2*abs(ewfd.Ey)^2 + ewfd.nzz^2*abs(ewfd.Ez)^2))), hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO*V_bias^2), (eps0*es.epsilonryy*wWG/hOEO*V_bias^2))');
model.result.numerical('int1').set('table', 'tbl7');
model.result.numerical('int1').selection.set([2]);
model.result.table.create('tbl8', 'Table');
model.result.table('tbl8').comments('Surface Integration 1 ((eps0*(ewfd.nxx^2*(abs(ewfd.Ex)^2 + ewfd.nyy^2*abs(ewfd.Ey)^2 + ewfd.nzz^2*abs(ewfd.Ez)^2))), hbar*omega*2/(eps0*es.epsilonryy*wWG/hOEO*V_bias^2), (eps0*es.epsilonryy*wWG/hOEO*V_bias^2))');
model.result.numerical('int1').set('table', 'tbl8');
model.result.numerical('int1').setResult;
model.result.numerical('int1').setIndex('expr', '(eps0*es.epsilonryy*wWG^2/hOEO*V_bias^2)', 2);
model.result.numerical('int1').setIndex('expr', '(eps0*es.epsilonryy*wWG^2/hOEO*V_bias^2)/[N]', 2);
model.result.numerical('int1').setIndex('expr', '', 2);
model.result.numerical('int1').setIndex('expr', '', 1);
model.result.numerical('int1').setIndex('expr', '(omega/4*eps0*(ewfd.nyy^2*((ewfd.Ey)*conj((ewfd.Ey))*es.Ey)', 0);
model.result.numerical('int1').setIndex('expr', '(omega/4*eps0*ewfd.nyy^2*ewfd.Ey*conj(ewfd.Ey)*es.Ey)', 0);
model.result.numerical('int1').setIndex('expr', '(omega/2*eps0*ewfd.nyy^4*ewfd.Ey*conj(ewfd.Ey)*es.Ey)*r', 0);
model.result.numerical('int1').setIndex('expr', '(omega/2*eps0*ewfd.nyy^4*ewfd.Ey*conj(ewfd.Ey)*es.Ey)*r33', 0);
model.result.numerical('int1').setIndex('expr', 'es.Ey', 1);
model.result.table.create('tbl9', 'Table');
model.result.table('tbl9').comments('Surface Integration 1 ((omega/2*eps0*ewfd.nyy^4*ewfd.Ey*conj(ewfd.Ey)*es.Ey)*r33, es.Ey, )');
model.result.numerical('int1').set('table', 'tbl9');
model.result.numerical('int1').setResult;

model.study('std1').feature('mode').setIndex('activate', false, 3);
model.study('std1').feature('mode').setIndex('activate', true, 3);
model.study('std1').feature('mode').setIndex('activate', false, 3);
model.study('std1').feature('mode').setIndex('activate', true, 3);
model.study('std1').feature('mode').setIndex('activate', false, 3);

model.sol('sol1').runAll;

model.result('pg1').run;
model.result('pg1').run;
model.result('pg1').feature('surf1').set('expr', 'es.Ey');
model.result('pg1').run;

model.study('std1').feature('mode').setIndex('activate', true, 3);

out = model;
