function out = model
%
% OEO_Slot.m
%
% Model exported on Feb 11 2020, 11:50 by COMSOL 5.4.0.246.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('Z:\Mitarbeiter\Haffner_Christian\oeffentlich\200124_Killian_Comsol');



% left over 
model.param.set('', 'sqrt(hbar)*2*pi*f_rf');

% Material properties
model.param.set('omega', '2*pi*c_const/wl');
model.param.set('gamma', '0.001*omega');
model.param.set('sigmaxx', '-1i*e^2/(4*pi*hbar)*log((-omega + 1i*2*gamma)/(+omega + 1i*2*gamma))');
model.param.set('sigmayy', '0 [S]');
model.param.set('sigmazz', 'sigmaxx');


model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').label('Thermal Oxide');
model.component('comp1').geom('geom1').feature('r1').set('pos', [0 0]);
model.component('comp1').geom('geom1').feature('r1').set('size', {'wSim' 'hSubstrate'});


model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').label('Waveguide');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'(wSim - wWG)/2' 'hSubstrate-hWG_top'});
model.component('comp1').geom('geom1').feature('r2').set('size', {'wWG' 'hWG_top'});
model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('r3').label('Organics');
model.component('comp1').geom('geom1').feature('r3').set('pos', {'0' 'hSubstrate'});
model.component('comp1').geom('geom1').feature('r3').set('size', {'wSim' 'hOrganic'});
model.component('comp1').geom('geom1').create('r4', 'Rectangle');
model.component('comp1').geom('geom1').feature('r4').label('GrapheneELectrodes');
model.component('comp1').geom('geom1').feature('r4').set('pos', {'0' 'hSubstrate-hWG/2-hContact/2'});
model.component('comp1').geom('geom1').feature('r4').set('size', {'wSim' 'hContact'});
model.component('comp1').geom('geom1').create('r5', 'Rectangle');
model.component('comp1').geom('geom1').feature('r5').label('Air');
model.component('comp1').geom('geom1').feature('r5').set('pos', {'0' 'hSubstrate+hOrganic+hBuffer'});
model.component('comp1').geom('geom1').feature('r5').set('size', {'wSim' '4e-6'});
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material.create('mat3', 'Common');
model.component('comp1').material.create('mat4', 'Common');
model.component('comp1').material.create('mat5', 'Common');
model.component('comp1').material.create('mat6', 'Common');
model.component('comp1').material.create('mat7', 'Common');
model.component('comp1').material('mat1').selection.set([2 7 9]);
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat2').selection.set([1 2 3 4 5 9 10]);
model.component('comp1').material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat3').selection.set([6 8]);
model.component('comp1').material('mat3').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat4').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat5').selection.set([7]);
model.component('comp1').material('mat5').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat6').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat7').selection.set([7]);
model.component('comp1').material('mat7').propertyGroup.create('RefractiveIndex', 'Refractive index');

model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
model.component('comp1').physics('ewfd').create('sctr1', 'Scattering', 1);
model.component('comp1').physics('ewfd').feature('sctr1').selection.set([1 2 3 4 5 6 7 8 9 10 11 15 17 21 23 24 25 26 27 28 29]);
model.component('comp1').physics('ewfd').create('scu1', 'SurfaceCurrent', 1);
model.component('comp1').physics('ewfd').feature('scu1').selection.set([6 15 17 21]);
model.component('comp1').physics.create('es', 'Electrostatics', 'geom1');
model.component('comp1').physics('es').selection.set([7]);
model.component('comp1').physics('es').create('gnd1', 'Ground', 1);
model.component('comp1').physics('es').feature('gnd1').selection.set([4 15 21]);
model.component('comp1').physics('es').create('pot1', 'ElectricPotential', 1);
model.component('comp1').physics('es').feature('pot1').selection.set([8 17 18 24]);

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri3', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri4', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri5', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri6', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').selection.set([2 7 9]);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([6 7 8]);
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').selection.set([6 8]);
model.component('comp1').mesh('mesh1').feature('ftri3').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri4').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri4').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri5').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri5').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').selection.set([4]);
model.component('comp1').mesh('mesh1').feature('ftri6').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri6').selection.set([1 2 3 4 5 9 10]);
model.component('comp1').mesh('mesh1').feature('ftri6').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').selection.set([1 3 10]);

model.result.table('tbl1').label('OEO Waveguide');
model.result.table('tbl1').comments('Global Evaluation 1 (10*log10(exp(4*pi/wl*imag(ewfd.neff)*1e-3)), real(ewfd.neff))');
model.result.table('tbl2').label('OEO Slot waveguide');
model.result.table('tbl2').comments('Global Evaluation 1 (10*log10(exp(4*pi/wl*imag(ewfd.neff)*1e-3)), real(ewfd.neff))');
model.result.table('tbl7').comments('g_0 (eps0*ewfd.epsilonrxx^2*r33*(ewfd.Ey*conj(ewfd.Ey))*c_const/wl*2*pi*V_bias/(ratio_wg*hWG))');

model.component('comp1').view('view1').axis.set('xmin', 1.1144627478643088E-6);
model.component('comp1').view('view1').axis.set('xmax', 5.291774868965149E-6);
model.component('comp1').view('view1').axis.set('ymin', 2.702558731471072E-6);
model.component('comp1').view('view1').axis.set('ymax', 4.948265996063128E-6);

model.component('comp1').material('mat1').label('Au');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'real(sqrt(1-wl^2/lp^2+1i*Gamma*wl^3/(2*pi*c_const*lp^2)))' '0' '0' '0' 'real(sqrt(1-wl^2/lp^2+1i*Gamma*wl^3/(2*pi*c_const*lp^2)))' '0' '0' '0' 'real(sqrt(1-wl^2/lp^2+1i*Gamma*wl^3/(2*pi*c_const*lp^2)))'});
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', {'imag(sqrt(1-wl^2/lp^2+1i*Gamma*wl^3/(2*pi*c_const*lp^2)))' '0' '0' '0' 'imag(sqrt(1-wl^2/lp^2+1i*Gamma*wl^3/(2*pi*c_const*lp^2)))' '0' '0' '0' 'imag(sqrt(1-wl^2/lp^2+1i*Gamma*wl^3/(2*pi*c_const*lp^2)))'});
model.component('comp1').material('mat2').label('SiO2');
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'1.44*1.44' '0' '0' '0' '1.44*1.44' '0' '0' '0' '1.44*1.44'});
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('n', {'sqrt(1+ASiO2*wl^2/(wl^2-l0SiO2^2) + BSiO2*wl^2/(wl^2-l1SiO2^2) + CSiO2*wl^2/(wl^2-l2SiO2^2))' '0' '0' '0' 'sqrt(1+ASiO2*wl^2/(wl^2-l0SiO2^2) + BSiO2*wl^2/(wl^2-l1SiO2^2) + CSiO2*wl^2/(wl^2-l2SiO2^2))' '0' '0' '0' 'sqrt(1+ASiO2*wl^2/(wl^2-l0SiO2^2) + BSiO2*wl^2/(wl^2-l1SiO2^2) + CSiO2*wl^2/(wl^2-l2SiO2^2))'});
model.component('comp1').material('mat2').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.component('comp1').material('mat3').label('Si');
model.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('n', {'real(sqrt(eps_Si_im(wl)+0*i*eps_Si_re(wl)))' '0' '0' '0' 'real(sqrt(eps_Si_im(wl)+0*i*eps_Si_re(wl)))' '0' '0' '0' 'real(sqrt(eps_Si_im(wl)+0*i*eps_Si_re(wl)))'});
model.component('comp1').material('mat3').propertyGroup('RefractiveIndex').set('ki', {'imag(sqrt(eps_Si_im(wl)+0*i*eps_Si_re(wl)))' '0' '0' '0' 'imag(sqrt(eps_Si_im(wl)+0*i*eps_Si_re(wl)))' '0' '0' '0' 'imag(sqrt(eps_Si_im(wl)+0*i*eps_Si_re(wl)))'});
model.component('comp1').material('mat4').label('SiNx');
model.component('comp1').material('mat4').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat4').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat4').propertyGroup('RefractiveIndex').set('n', {'sqrt(1+ASiNx*wl^2/(wl^2-l0SiNx^2) + BSiNx*wl^2/(wl^2-l1SiNx^2))' '0' '0' '0' 'sqrt(1+ASiNx*wl^2/(wl^2-l0SiNx^2) + BSiNx*wl^2/(wl^2-l1SiNx^2))' '0' '0' '0' 'sqrt(1+ASiNx*wl^2/(wl^2-l0SiNx^2) + BSiNx*wl^2/(wl^2-l1SiNx^2))'});
model.component('comp1').material('mat4').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.component('comp1').material('mat5').label('Organics');
model.component('comp1').material('mat5').propertyGroup('def').set('relpermittivity', {'3' '0' '0' '0' '3' '0' '0' '0' '3'});
model.component('comp1').material('mat5').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat5').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat5').propertyGroup('RefractiveIndex').set('n', {'real(sqrt(epsrOrganics(wl)+0*i*epsiOrganics(wl)))' '0' '0' '0' 'real(sqrt(epsrOrganics(wl)+0*i*epsiOrganics(wl)))' '0' '0' '0' 'real(sqrt(epsrOrganics(wl)+0*i*epsiOrganics(wl)))'});
model.component('comp1').material('mat5').propertyGroup('RefractiveIndex').set('ki', {'imag(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))' '0' '0' '0' 'imag(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))' '0' '0' '0' 'imag(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))'});
model.component('comp1').material('mat6').label('Al2O3');
model.component('comp1').material('mat6').propertyGroup('def').set('relpermittivity', {'9.8' '0' '0' '0' '9.8' '0' '0' '0' '9.8'});
model.component('comp1').material('mat6').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat6').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat6').propertyGroup('RefractiveIndex').set('n', {'sqrt(1+AAl2O3*wl^2/(wl^2-l0Al2O3^2) + BAl2O3*wl^2/(wl^2-l1Al2O3^2) + CAl2O3*wl^2/(wl^2-l2Al2O3^2))' '0' '0' '0' 'sqrt(1+AAl2O3*wl^2/(wl^2-l0Al2O3^2) + BAl2O3*wl^2/(wl^2-l1Al2O3^2) + CAl2O3*wl^2/(wl^2-l2Al2O3^2))' '0' '0' '0' 'sqrt(1+AAl2O3*wl^2/(wl^2-l0Al2O3^2) + BAl2O3*wl^2/(wl^2-l1Al2O3^2) + CAl2O3*wl^2/(wl^2-l2Al2O3^2))'});
model.component('comp1').material('mat6').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.component('comp1').material('mat7').label('Organics ES');
model.component('comp1').material('mat7').propertyGroup('def').set('relpermittivity', {'3' '0' '0' '0' '3' '0' '0' '0' '3'});
model.component('comp1').material('mat7').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat7').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat7').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat7').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat7').propertyGroup('RefractiveIndex').set('n', {'real(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))' '0' '0' '0' 'real(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))-0.5*real(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))^3*r33*es.Ey' '0' '0' '0' 'real(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))'});
model.component('comp1').material('mat7').propertyGroup('RefractiveIndex').set('ki', {'imag(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))' '0' '0' '0' 'imag(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))' '0' '0' '0' 'imag(sqrt(epsrOrganics(wl)+i*epsiOrganics(wl)))'});

model.component('comp1').physics('ewfd').feature('scu1').set('Js0', {'sigmaxx*ewfd.Ex'; 'sigmayy*ewfd.Ey'; 'sigmazz*ewfd.Ez'});
model.component('comp1').physics('es').feature('pot1').set('V0', 'V_bias');

model.component('comp1').mesh('mesh1').feature('ftri1').set('yscale', 5);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmax', '50[nm]');
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmax', 'wl/(3.4*20)');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hmax', 'wl/(2.0*20)');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmax', 'wl/(1.78*20)');
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftri4').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hmax', 'wl/(1.78*20)[nm]');
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftri5').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hmax', 'wl/(1.45*10)');
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hmin', 3.86E-11);
model.component('comp1').mesh('mesh1').feature('ftri6').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.study('std1').create('mode', 'ModeAnalysis');
model.study('std1').feature('stat').set('activate', {'ewfd' 'off' 'es' 'on' 'frame:spatial1' 'on'});

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').create('su1', 'StoreSolution');
model.sol('sol1').create('st2', 'StudyStep');
model.sol('sol1').create('v2', 'Variables');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol.create('sol3');
model.sol('sol3').study('std1');
model.sol('sol3').label('Parametric Solutions 1');

model.result.numerical.create('int1', 'IntSurface');
model.result.numerical.create('int2', 'IntSurface');
model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical('int1').set('data', 'dset3');
model.result.numerical('int1').selection.set([7]);
model.result.numerical('int1').set('probetag', 'none');
model.result.numerical('int2').set('data', 'dset3');
model.result.numerical('int2').selection.set([1 2 3 4 5 6 7 8 9 10]);
model.result.numerical('int2').set('probetag', 'none');
model.result.numerical('gev1').set('data', 'dset3');
model.result.numerical('gev1').set('probetag', 'none');
model.result.create('pg5', 'PlotGroup1D');
model.result.create('pg6', 'PlotGroup1D');
model.result.create('pg9', 'PlotGroup2D');
model.result.create('pg10', 'PlotGroup2D');
model.result.create('pg11', 'PlotGroup2D');
model.result.create('pg12', 'PlotGroup2D');
model.result('pg5').create('tblp1', 'Table');
model.result('pg6').create('tblp1', 'Table');
model.result('pg9').set('data', 'dset3');
model.result('pg9').create('surf1', 'Surface');
model.result('pg10').set('data', 'dset3');
model.result('pg10').create('surf1', 'Surface');
model.result('pg11').create('surf1', 'Surface');
model.result('pg12').create('surf1', 'Surface');

model.study('std1').feature('mode').set('modeFreq', 'c_const/wl');
model.study('std1').feature('mode').set('neigsactive', true);
model.study('std1').feature('mode').set('neigs', 2);
model.study('std1').feature('mode').set('shiftactive', true);
model.study('std1').feature('mode').set('shift', '3.4');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st2').set('studystep', 'mode');
model.sol('sol1').feature('v2').set('initmethod', 'sol');
model.sol('sol1').feature('v2').set('initsol', 'sol1');
model.sol('sol1').feature('v2').set('initsoluse', 'sol2');
model.sol('sol1').feature('v2').set('solnum', 'auto');
model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
model.sol('sol1').feature('v2').set('notsol', 'sol1');
model.sol('sol1').feature('v2').set('notsolnum', 'auto');
model.sol('sol1').feature('e1').set('control', 'mode');
model.sol('sol1').feature('e1').set('transform', 'effective_mode_index');
model.sol('sol1').feature('e1').set('neigs', 2);
model.sol('sol1').feature('e1').set('shift', '3.4');
model.sol('sol1').feature('e1').set('linpmethod', 'sol');
model.sol('sol1').feature('e1').set('linpsol', 'sol1');
model.sol('sol1').feature('e1').set('linpsoluse', 'sol2');
model.sol('sol1').feature('e1').set('solnum', 'auto');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').runAll;

model.result.numerical('int1').label('g_0');
model.result.numerical('int1').set('table', 'tbl7');
model.result.numerical('int1').set('expr', {'eps0*ewfd.epsilonrxx^2*r33*(ewfd.Ey*conj(ewfd.Ey))*c_const/wl*2*pi*V_bias/(ratio_wg*hWG)'});
model.result.numerical('int1').set('unit', {'N/s'});
model.result.numerical('int1').set('descr', {'Vaccum Coupling'});
model.result.numerical('int2').label('norm Field');
model.result.numerical('int2').set('table', 'tbl7');
model.result.numerical('int2').set('expr', {'V_bias*2*eps0*ewfd.epsilonrxx*((ewfd.Ex*conj(ewfd.Ex)) + (ewfd.Ey*conj(ewfd.Ey)) + (ewfd.Ez*conj(ewfd.Ez)))'});
model.result.numerical('int2').set('unit', {'kg^2*m^3/(s^5*A)'});
model.result.numerical('int2').set('descr', {'norm Field'});
model.result.numerical('gev1').set('table', 'tbl7');
model.result.numerical('gev1').set('expr', {'ewfd.neff'});
model.result.numerical('gev1').set('unit', {'1'});
model.result.numerical('gev1').set('descr', {'Effective mode index'});
model.result.numerical('int1').setResult;
model.result.numerical('int2').appendResult;
model.result.numerical('gev1').appendResult;
model.result('pg5').label('OEO Waveguide');
model.result('pg5').set('data', 'none');
model.result('pg5').set('xlabel', 'ratio_wg');
model.result('pg5').set('xlabelactive', true);
model.result('pg5').set('ylabel', 'real(ewfd.neff) (1)');
model.result('pg5').set('ylabelactive', true);
model.result('pg5').set('axislimits', true);
model.result('pg5').set('xmin', 0.03418930549072833);
model.result('pg5').set('xmax', 0.9584585527695538);
model.result('pg5').set('ymin', 1.5445);
model.result('pg5').set('ymax', 1.5454);
model.result('pg5').set('showlegends', false);
model.result('pg5').feature('tblp1').set('xaxisdata', 2);
model.result('pg5').feature('tblp1').set('plotcolumninput', 'manual');
model.result('pg5').feature('tblp1').set('linestyle', 'none');
model.result('pg5').feature('tblp1').set('linemarker', 'cyclereset');
model.result('pg5').feature('tblp1').set('markerpos', 'datapoints');
model.result('pg5').feature('tblp1').set('legend', true);
model.result('pg5').feature('tblp1').set('legendmethod', 'manual');
model.result('pg5').feature('tblp1').set('legends', {'db/mm'});
model.result('pg6').label('OEO Slot waveguide');
model.result('pg6').set('data', 'none');
model.result('pg6').set('xlabel', 'ratio_wg');
model.result('pg6').set('ylabel', 'real(ewfd.neff) (1)');
model.result('pg6').set('xlabelactive', false);
model.result('pg6').set('ylabelactive', false);
model.result('pg6').feature('tblp1').set('table', 'tbl2');
model.result('pg6').feature('tblp1').set('xaxisdata', 2);
model.result('pg6').feature('tblp1').set('plotcolumninput', 'manual');
model.result('pg6').feature('tblp1').set('linestyle', 'none');
model.result('pg6').feature('tblp1').set('linemarker', 'cycle');
model.result('pg6').feature('tblp1').set('markerpos', 'datapoints');
model.result('pg9').label('Electric Field (ewfd) 1');
model.result('pg9').set('looplevel', [2 1 2]);
model.result('pg9').set('frametype', 'spatial');
model.result('pg9').feature('surf1').set('expr', 'abs(ewfd.Ey)');
model.result('pg9').feature('surf1').set('descr', 'abs(ewfd.Ey)');
model.result('pg9').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg9').feature('surf1').set('smooth', 'internal');
model.result('pg9').feature('surf1').set('resolution', 'normal');
model.result('pg10').label('Electric Potential (es) 1');
model.result('pg10').set('looplevel', [1 2 1]);
model.result('pg10').set('frametype', 'spatial');
model.result('pg10').feature('surf1').set('expr', 'abs(es.Ey)');
model.result('pg10').feature('surf1').set('descr', 'abs(es.Ey)');
model.result('pg10').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg10').feature('surf1').set('resolution', 'normal');
model.result('pg11').label('Electric Field (ewfd)');
model.result('pg11').set('frametype', 'spatial');
model.result('pg11').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg11').feature('surf1').set('smooth', 'internal');
model.result('pg11').feature('surf1').set('resolution', 'normal');
model.result('pg12').label('Electric Potential (es)');
model.result('pg12').set('frametype', 'spatial');
model.result('pg12').feature('surf1').set('expr', 'V');
model.result('pg12').feature('surf1').set('unit', 'V');
model.result('pg12').feature('surf1').set('descr', 'Electric potential');
model.result('pg12').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg12').feature('surf1').set('resolution', 'normal');


model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl7', 'Table');


out = model;
