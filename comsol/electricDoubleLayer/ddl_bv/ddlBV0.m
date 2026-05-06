function out = model
%
% ddlBV.m
%
% Model exported on Oct 10 2025, 11:56 by COMSOL 6.2.0.290.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/howardtu/Documents/modeling/continuum/feecm/comsol/ddl_bv/ref');

model.label('diffuse_double_layer_with_charge_transfer.mph');

model.title('Diffuse Double Layer with Charge Transfer');

model.description(['At the electrode-electrolyte interface, there is a thin layer of space charge, called the diffuse double layer. In this region, electroneutrality does not hold. The double layer may be of interest when modeling devices such as electrochemical supercapacitors and nanoelectrodes.' newline  newline 'The Diffuse Double Layer tutorial shows you how to couple the Nernst' native2unicode(hex2dec({'20' '13'}), 'unicode') 'Planck equations to the Poisson equation in order to describe a diffuse double layer according to the Gouy' native2unicode(hex2dec({'20' '13'}), 'unicode') 'Chapman' native2unicode(hex2dec({'20' '13'}), 'unicode') 'Stern model.' newline  newline 'The simulation application extends the simple example by including two electrodes. It also considers Faradaic (charge transfer) electrode reactions. An additional equation is solved to ensure overall conservation of charge.']);

model.param.set('epsilon', '0.001', 'Dimensionless Debye length scale');
model.param.set('L', 'lambdaD/epsilon', 'Cell length');
model.param.set('Dp', '1e-9[m^2/s]', 'Diffusion coefficient, positive ion');
model.param.set('T', '298.15[K]', 'Temperature');
model.param.set('RT', 'R_const*T', 'Molar gas constant * Temperature');
model.param.set('Dm', 'Dp', 'Diffusion coefficient, negative ion');
model.param.set('cref', '1[mol/m^3]', 'Reference ion concentration');
model.param.set('cM', '1[mol/m^3]', 'Metal reference concentration');
model.param.set('Z_ch', '1', 'Ion charge');
model.param.set('alphac', '0.5', 'Cathodic charge transfer coefficient');
model.param.set('alphaa', '1-alphac', 'Anodic charge transfer coefficient');
model.param.set('jr', '10', 'Dimensionless anodic reaction current density');
model.param.set('kc', '10', 'Dimensionless cathodic rate coefficient');
model.param.set('delta', '0.025', 'Dimensionless Stern layer thickness');
model.param.set('J', '0.9', 'Dimensionless cell current density');
model.param.set('Kc', 'kc*4*Dp/L', 'Cathodic rate constant');
model.param.set('Ka', 'jr*4*Dp*cref/(L*cM)', 'Anodic rate constant');
model.param.set('id', '4*Z_ch*F_const*Dp*cref/L', 'Nernst limiting current density');
model.param.set('icell', 'J*id', 'Cell current density');
model.param.set('lambdaD', 'sqrt(epsilon0_const*eps_H2O*RT/(2*Z_ch^2*F_const^2*cref))', 'Debye length');
model.param.set('lambdaS', 'delta*lambdaD', 'Stern layer thickness');
model.param.set('eps_H2O', '78.5', 'Relative permittivity of water');
model.param.set('eps_S', '10', 'Relative permittivity of Stern layer');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 1);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').set('coord', {'0' 'L'});
model.component('comp1').geom('geom1').run;

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('deltaphi', 'tcd.phisext-phil', 'Metal - reaction plane potential difference');
model.component('comp1').variable('var1').set('rho_s', 'epsilon0_const*eps_S*deltaphi/lambdaS', 'Surface charge density');

model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').selection.set([1]);

model.component('comp1').physics.create('tcd', 'TertiaryCurrentDistributionNernstPlanck', 'geom1');
model.component('comp1').physics('tcd').field('concentration').field('cp');
model.component('comp1').physics('tcd').field('concentration').component({'cp' 'cm'});
model.component('comp1').physics('tcd').create('es1', 'ElectrodeSurface', 0);
model.component('comp1').physics('tcd').feature('es1').selection.set([1]);
model.component('comp1').physics('tcd').create('es2', 'ElectrodeSurface', 0);
model.component('comp1').physics('tcd').feature('es2').selection.set([2]);
model.component('comp1').physics('tcd').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('tcd').feature('sfcd1').selection.all;
model.component('comp1').physics('tcd').create('gconstr1', 'GlobalConstraint', -1);

model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').create('size2', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').selection.geom('geom1', 0);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').selection.all;

model.component('comp1').view('view1').axis.set('xmin', -4.809917300008237E-7);
model.component('comp1').view('view1').axis.set('xmax', 1.0100819963554386E-5);

model.common('cminpt').set('modified', {'temperature' 'T'});

model.component('comp1').physics('tcd').prop('SpeciesProperties').set('ChargeTransportModel', 'Poisson');
model.component('comp1').physics('tcd').feature('sp1').set('z', {'Z_ch'; '-Z_ch'});
model.component('comp1').physics('tcd').feature('ice1').set('D_cp', {'Dp'; '0'; '0'; '0'; 'Dp'; '0'; '0'; '0'; 'Dp'});
model.component('comp1').physics('tcd').feature('ice1').set('D_cm', {'Dm'; '0'; '0'; '0'; 'Dm'; '0'; '0'; '0'; 'Dm'});
model.component('comp1').physics('tcd').feature('init1').set('initc', {'cref'; 'cref'});
model.component('comp1').physics('tcd').feature('es1').feature('er1').set('Vi0', [-1; 0]);
model.component('comp1').physics('tcd').feature('es1').feature('er1').set('Eeq_mat', 'userdef');
model.component('comp1').physics('tcd').feature('es1').feature('er1').set('ElectrodeKinetics', 'ConcentrationDependentKinetics');
model.component('comp1').physics('tcd').feature('es1').feature('er1').set('i0', 'Ka*F_const*cM*Z_ch');
model.component('comp1').physics('tcd').feature('es1').feature('er1').set('alphaa', 'alphaa*Z_ch');
model.component('comp1').physics('tcd').feature('es1').feature('er1').set('alphac', 'alphac*Z_ch');
model.component('comp1').physics('tcd').feature('es1').feature('er1').set('CO', 'Kc/Ka*cp/cM');
model.component('comp1').physics('tcd').feature('es2').set('BoundaryCondition', 'AverageCurrentDensity');
model.component('comp1').physics('tcd').feature('es2').set('Ial', 'icell');
model.component('comp1').physics('tcd').feature('es2').feature('er1').set('Vi0', [-1; 0]);
model.component('comp1').physics('tcd').feature('es2').feature('er1').set('Eeq_mat', 'userdef');
model.component('comp1').physics('tcd').feature('es2').feature('er1').set('ElectrodeKinetics', 'ConcentrationDependentKinetics');
model.component('comp1').physics('tcd').feature('es2').feature('er1').set('i0', 'Ka*F_const*cM*Z_ch');
model.component('comp1').physics('tcd').feature('es2').feature('er1').set('alphaa', 'alphaa*Z_ch');
model.component('comp1').physics('tcd').feature('es2').feature('er1').set('alphac', 'alphac*Z_ch');
model.component('comp1').physics('tcd').feature('es2').feature('er1').set('CO', 'Kc/Ka*cp/cM');
model.component('comp1').physics('tcd').feature('sfcd1').set('rhos', 'rho_s');
model.component('comp1').physics('tcd').feature('gconstr1').set('constraintExpression', 'intop1(cm)-(cref*L)');

model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 'L/20');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmin', 2.89E-9);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hmax', 'lambdaD/10');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hmin', 2.89E-9);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hminactive', false);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('param', 'Parametric');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').create('i2', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol.create('sol2');
model.sol('sol2').study('std1');
model.sol('sol2').label('Parametric Solutions 1');

model.batch.create('p1', 'Parametric');
model.batch('p1').create('so1', 'Solutionseq');
model.batch('p1').study('std1');

model.result.create('pg1', 'PlotGroup1D');
model.result.create('pg2', 'PlotGroup1D');
model.result.create('pg3', 'PlotGroup1D');
model.result('pg1').set('data', 'dset2');
model.result('pg1').create('lngr1', 'LineGraph');
model.result('pg1').feature('lngr1').set('xdata', 'expr');
model.result('pg1').feature('lngr1').selection.set([1]);
model.result('pg1').feature('lngr1').set('expr', '(cp+cm)/(2*cref)');
model.result('pg2').set('data', 'dset2');
model.result('pg2').create('lngr1', 'LineGraph');
model.result('pg2').feature('lngr1').set('xdata', 'expr');
model.result('pg2').feature('lngr1').selection.set([1]);
model.result('pg2').feature('lngr1').set('expr', '(cp-cm)/(2*cref)');
model.result('pg3').set('data', 'dset2');
model.result('pg3').create('lngr1', 'LineGraph');
model.result('pg3').feature('lngr1').set('xdata', 'expr');
model.result('pg3').feature('lngr1').selection.set([1]);
model.result('pg3').feature('lngr1').set('expr', 'phil*Z_ch*F_const/(R_const*T)');

model.study('std1').feature('param').set('pname', {'epsilon'});
model.study('std1').feature('param').set('plistarr', {'0.001 0.01 0.1'});
model.study('std1').feature('param').set('punit', {''});

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Stationary');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').feature('comp1_phil').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_phil').set('scaleval', 1);
model.sol('sol1').feature('v1').feature('comp1_tcd_es2_phisext').set('scalemethod', 'manual');
model.sol('sol1').feature('v1').feature('comp1_tcd_es2_phisext').set('scaleval', 1);
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').set('stol', 1.0E-4);
model.sol('sol1').feature('s1').set('probesel', 'none');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'd1');
model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('s1').feature('d1').label('Direct (tcd)');
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i1').label('Algebraic Multigrid (tcd)');
model.sol('sol1').feature('s1').feature('i1').set('nlinnormuse', true);
model.sol('sol1').feature('s1').feature('i1').set('maxlinit', 1000);
model.sol('sol1').feature('s1').feature('i1').feature('ilDef').label('Incomplete LU 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').label('Multigrid 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('compactaggregation', true);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i2').label('Geometric Multigrid (tcd)');
model.sol('sol1').feature('s1').feature('i2').set('nlinnormuse', true);
model.sol('sol1').feature('s1').feature('i2').set('maxlinit', 1000);
model.sol('sol1').feature('s1').feature('i2').feature('ilDef').label('Incomplete LU 1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').label('Multigrid 1.1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('pr').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('po').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol1').feature('s1').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').runAll;

model.batch('p1').set('pname', {'epsilon'});
model.batch('p1').set('plistarr', {'0.001 0.01 0.1'});
model.batch('p1').set('punit', {''});
model.batch('p1').set('err', true);
model.batch('p1').feature('so1').set('seq', 'sol1');
model.batch('p1').feature('so1').set('store', true);
model.batch('p1').feature('so1').set('psol', 'sol2');
model.batch('p1').feature('so1').set('keeprom', true);
model.batch('p1').feature('so1').set('param', {'"epsilon","0.001"' '"epsilon","0.01"' '"epsilon","0.1"'});
model.batch('p1').attach('std1');
model.batch('p1').run;

model.result('pg1').label('Dimensionless concentration');
model.result('pg1').set('titletype', 'label');
model.result('pg1').set('xlabel', 'x/L (1)');
model.result('pg1').set('ylabel', '(cp+cm)/(2*cref) (1)');
model.result('pg1').set('legendpos', 'lowerright');
model.result('pg1').set('xlabelactive', false);
model.result('pg1').set('ylabelactive', false);
model.result('pg1').feature('lngr1').set('xdataexpr', 'x/L');
model.result('pg1').feature('lngr1').set('xdataunit', '1');
model.result('pg1').feature('lngr1').set('xdatadescr', 'x/L');
model.result('pg1').feature('lngr1').set('linewidth', 'preference');
model.result('pg1').feature('lngr1').set('legend', true);
model.result('pg1').feature('lngr1').set('resolution', 'normal');
model.result('pg2').label('Dimensionless charge density');
model.result('pg2').set('titletype', 'label');
model.result('pg2').set('xlabel', 'x/L (1)');
model.result('pg2').set('ylabel', '(cp-cm)/(2*cref) (1)');
model.result('pg2').set('xlabelactive', false);
model.result('pg2').set('ylabelactive', false);
model.result('pg2').feature('lngr1').set('xdataexpr', 'x/L');
model.result('pg2').feature('lngr1').set('xdataunit', '1');
model.result('pg2').feature('lngr1').set('xdatadescr', 'x/L');
model.result('pg2').feature('lngr1').set('linewidth', 'preference');
model.result('pg2').feature('lngr1').set('legend', true);
model.result('pg2').feature('lngr1').set('resolution', 'normal');
model.result('pg3').label('Dimensionless potential');
model.result('pg3').set('titletype', 'label');
model.result('pg3').set('xlabel', 'x/L (1)');
model.result('pg3').set('ylabel', 'phil*Z_ch*F_const/(R_const*T) (1)');
model.result('pg3').set('legendpos', 'lowerright');
model.result('pg3').set('xlabelactive', false);
model.result('pg3').set('ylabelactive', false);
model.result('pg3').feature('lngr1').set('xdataexpr', 'x/L');
model.result('pg3').feature('lngr1').set('xdataunit', '1');
model.result('pg3').feature('lngr1').set('xdatadescr', 'x/L');
model.result('pg3').feature('lngr1').set('linewidth', 'preference');
model.result('pg3').feature('lngr1').set('legend', true);
model.result('pg3').feature('lngr1').set('resolution', 'normal');

out = model;
