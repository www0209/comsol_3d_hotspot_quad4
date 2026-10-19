function out = t3d_run(varargin)
% T3D_RUN  (comsol_3d_hotspot)  Three-dimensional COMSOL model of a StackEM stack: stress on every rail with the
% HOTSPOT temperatures.  Same geometry, rails, equations, outputs and pictures as comsol_3d\t3d_run.m; the only
% difference is that the temperature is NOT solved in COMSOL: the rails take T_L, T_R of StackEM's tables (the HotSpot
% field) and the temperature pictures show the HotSpot field painted on the dies.  There is no thermal study, so no
% comparison of COMSOL against HotSpot temperatures: the two temperature columns of the csv files are both HotSpot.
% The case (base3 / quad4 / dense3 / ev6) is whatever data\<case>\stack.json describes; the script itself is generic.
%
%   t3d_run('max_rails', 2, 'tag', 'smoke')      first test
%   t3d_run                                      all rails, 10 years
%
% OPTIONS (name, value)                        default
%   'tag'         output folder out_<tag>          'run'
%   'tend_years'  end of the stress integration    10
%   'dt_decade'   output spacing in log10(t)       0.05
%   'horizons'    years at which fields are saved  [1 2 5 10]
%   'max_rails'   per die and direction            inf
%   'nper'        line elements per segment        60
%   'order'       element order, stress            1
%   'hdie', 'hspr', 'hsink', 'nz', 'zscale'      volume mesh of the package (pictures only here; no equation is solved
%                                                 in the solid)   stack.json mesh.pictures  (coarse)
%   'inserts'     1 = draw the inserts of stack.json as geometry (ev6: TSV strips)   1
%   'rtol'        relative tolerance              1e-6
%   'field_pts'   points per segment in t3d_field_<H>yr.csv   11
%   'rail_plots'  one 1-D stress plot per rail exported by COMSOL   1
%   'topview'     per-die pictures from a camera above the die      1
%   'history'     1 = save t3d_history.mat (time x position field; ev6: about 2 GB)   1
%   'data'        <this folder>\data\<case>        the only data folder next to this file
%
% OUTPUT (out_<tag>\): the same files as comsol_3d (<case>_3d.mph, t3d_die_temperature.csv, t3d_nodes.csv,
% t3d_nodes_history.csv, t3d_history.mat, t3d_field_<H>yr.csv, t3d_rail_result.csv, t3d_rail_smax.csv, t3d_meta.json,
% comsol_images\).  t3d_figures.py and t3d_compare.py work unchanged.

import com.comsol.model.*
import com.comsol.model.util.*

here = fileparts(mfilename('fullpath'));  Y = 365.25 * 86400;
dd_ = dir(fullfile(here, 'data'));  dd_ = dd_([dd_.isdir] & ~ismember({dd_.name}, {'.', '..'}));
if numel(dd_) ~= 1, error('t3d:data', 'expected exactly one case folder under %s, found %d', fullfile(here, 'data'), numel(dd_)); end
o = struct('tag','run', 'tsource','hotspot', 'tend_years',10, 'dt_decade',0.05, 'horizons',[1 2 5 10], 'max_rails',inf, 'nper',60, ...
           'order',1, 'order_T',1, 'hdie',[], 'hspr',[], 'hsink',[], 'nz',[], 'zscale','true', 'mesh3d',4, 'inserts',1, ...
           'rtol',1e-6, 'field_pts',11, 'rail_plots',1, 'topview',1, 'history',1, 'data', fullfile(here, 'data', dd_(1).name));
for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end
useHS = true;  truez = strcmpi(o.zscale, 'true');
outDir = fullfile(here, ['out_' o.tag]);  if ~exist(outDir, 'dir'), mkdir(outDir); end
imgDir = fullfile(outDir, 'comsol_images');  if ~exist(imgDir, 'dir'), mkdir(imgDir); end
fb = {};

R  = dlmread(fullfile(o.data, 'rails.csv'), ',', 1, 0);       % rail,die_id,kind_is_V,index,n_seg,x0,y0,W,H,Gamma,immortal,n_cells,sigma_ss,seg_offset
S  = dlmread(fullfile(o.data, 'segments.csv'), ',', 1, 0);    % rail,k,L,j,T_L,T_R,T_m,Gamma,sigma_T
PB = dlmread(fullfile(o.data, 'power_blocks.csv'), ',', 1, 0);% die_id,x0,x1,y0,y1,P
HT = dlmread(fullfile(o.data, 'hotspot_T.csv'), ',', 1, 0);   % die_id,x,y,T
C  = jsondecode(fileread(fullfile(o.data, 'constants.json')));
K  = jsondecode(fileread(fullfile(o.data, 'stack.json')));
caseName = K.case;  dieNames = K.die_names;  if ~iscell(dieNames), dieNames = cellstr(dieNames); end
mesh_def = K.mesh.pictures;                                    % coarse: nothing is solved in the solid (see t3d_make_data.py)
if isempty(o.hdie), o.hdie = mesh_def.hdie; end;  if isempty(o.hspr), o.hspr = mesh_def.hspr; end
if isempty(o.hsink), o.hsink = mesh_def.hsink; end;  if isempty(o.nz), o.nz = mesh_def.nz(:)'; end
% ---------------- output times and horizons ----------------
hz = sort(unique([o.horizons(o.horizons <= o.tend_years), o.tend_years]));  nh = numel(hz);
tl = unique(round([0, 10 .^ (4:o.dt_decade:log10(o.tend_years * Y)), hz * Y], 3));  nt = numel(tl);
ih = zeros(1, nh);  for h = 1:nh, [~, ih(h)] = min(abs(tl - hz(h) * Y)); end
sc = C.em.sigma_crit;

% ---------------- rails kept ----------------
keep = false(size(R, 1), 1);
for d = unique(R(:, 2))'
    for kv = 0:1
        m = find(R(:, 2) == d & R(:, 3) == kv);  [~, s] = sort(R(m, 4));  m = m(s);
        keep(m(1:min(numel(m), o.max_rails))) = true;
    end
end
idx = find(keep);  nr = numel(idx);  Sr = R(idx, 5);  W = R(idx, 8);  die = R(idx, 2);  isV = R(idx, 3) == 1;
x0r = R(idx, 6);  y0r = R(idx, 7);  nd = numel(dieNames);

% ---------------- layers of the stack (bottom -> top, from stack.json) ----------------
tm = 50e-6;                                                   % drawn thickness of bond / TIM in the 'scaled' model
LJ = K.layers;  if isstruct(LJ), LJ = num2cell(LJ); end       % jsondecode: struct array or cell array of structs
IJ = {};  if isfield(K, 'inserts') && ~isempty(K.inserts), IJ = K.inserts;  if isstruct(IJ), IJ = num2cell(IJ); end, end
L = struct('name', {}, 'kind', {}, 't', {}, 'tm', {}, 'k', {}, 'sx', {}, 'sy', {}, 'h', {}, 'nz', {}, 'die', {}, 'ins', {});
for i = 1:numel(LJ)
    j = LJ{i};  kind = lower(j.kind);  di = -1;
    switch kind
        case 'die',      h_ = o.hdie;   nz_ = o.nz(1);  di = j.die_index;  tm_ = j.thickness;
        case 'bond',     h_ = o.hdie;   nz_ = o.nz(2);  tm_ = tm;
        case 'tim',      h_ = o.hdie;   nz_ = o.nz(3);  tm_ = tm;
        case 'spreader', h_ = o.hspr;   nz_ = o.nz(4);  tm_ = j.thickness;
        case 'sink',     h_ = o.hsink;  nz_ = o.nz(5);  tm_ = j.thickness;
        otherwise, error('t3d:layer', 'unknown layer kind %s', kind);
    end
    ins = zeros(0, 5);
    for q = 1:numel(IJ)
        if strcmp(IJ{q}.layer, j.name), ins(end+1, :) = [IJ{q}.x0, IJ{q}.x1, IJ{q}.y0, IJ{q}.y1, IJ{q}.k]; end %#ok<AGROW>
    end
    L(end+1) = struct('name', j.name, 'kind', kind, 't', j.thickness, 'tm', tm_, 'k', j.k, 'sx', j.size_x, 'sy', j.size_y, 'h', h_, 'nz', nz_, 'die', di, 'ins', ins); %#ok<AGROW>
end
if truez, for i = 1:numel(L), L(i).tm = L(i).t; end, end
nl = numel(L);  zb = [0, cumsum([L.tm])];                     % zb(i) .. zb(i+1) = layer i
dieL = zeros(1, nd);                                          % layer index of every die (die_index 0 .. nd-1 -> dieL(1 .. nd))
for i = 1:nl, if L(i).die >= 0, dieL(L(i).die + 1) = i; end, end
if any(dieL == 0) || ~strcmp(L(nl).kind, 'sink') || ~strcmp(L(nl-1).kind, 'spreader'), error('t3d:layer', 'stack.json: every die needs a layer, the last two layers must be spreader and sink'); end
for d = 1:nd, if ~strcmp(L(dieL(d)).name, dieNames{d}), error('t3d:layer', 'die %d: layer %s but die name %s', d, L(dieL(d)).name, dieNames{d}); end, end
if any(diff(dieL) <= 0), error('t3d:layer', 'stack.json: die_index must increase bottom -> top (t3d_make_data.py numbers the dies that way)'); end
cd_ = C.dies;  if ~iscell(cd_), cd_ = cellstr(cd_); end
if ~isequal(cd_(:), dieNames(:)), error('t3d:layer', 'constants.json dies and stack.json die_names differ'); end
ztop = zb(dieL + 1);  tdie = [L(dieL).t];  zmid = ztop - 0.5 * tdie;    % top face and mid-plane of every die
dz = 2e-6;  zH = ztop + dz;  zV = ztop + 2 * dz;              % the two rail levels of every die (inside the layer above the die)
for d = 1:nd, if L(dieL(d) + 1).tm < 2.5 * dz, error('t3d:layer', 'layer %s above die %s is too thin for the rail levels', L(dieL(d) + 1).name, dieNames{d}); end, end
zr = zH(die + 1)';  zr(isV) = zV(die(isV) + 1)';
SX = K.die_size_x;  SY = K.die_size_y;  cx = 0.5 * SX;  cy = 0.5 * SY;   % all layers are centred on the die stack
sinkDoms = nl;  stackL = find(~ismember({L.kind}, {'spreader', 'sink'}));

% ---------------- tables: per-segment numbers against Sg = slot*Sstep + coordinate along the rail ----------------
t0 = tic;
aj = cell(nr, 1);  Lr = zeros(nr, 1);
for q = 1:nr
    Ls = S(R(idx(q), 14) + (1:Sr(q)), 3);  Lr(q) = sum(Ls);
    if isV(q), aj{q} = y0r(q) + [0; cumsum(Ls)]; else, aj{q} = x0r(q) + [0; cumsum(Ls)]; end
end
across = x0r .* isV + y0r .* ~isV;  slot = (0:nr-1)';
amax = max(cellfun(@(v) v(end), aj));  Sstep = 1.25 * amax + 10 * max(W);
ep = max(2e-9, 2e-11 * nr * Sstep);
tab = zeros(2 * sum(Sr), 9);  p = 0;                          % Sg, start of the segment, L, T_L, dT, T_m, lam, j, sigma_T
for q = 1:nr
    rows = R(idx(q), 14) + (1:Sr(q));  Ls = S(rows, 3);  a0 = aj{q};
    v = [a0(1:end-1), Ls, S(rows, 5), S(rows, 6) - S(rows, 5), S(rows, 7), Ls ./ S(rows, 8), S(rows, 4), S(rows, 9)];
    k = p + (1:2:2*Sr(q));
    tab(k, :) = [slot(q) * Sstep + a0(1:end-1) + ep, v];  tab(k + 1, :) = [slot(q) * Sstep + a0(2:end) - ep, v];
    p = p + 2 * Sr(q);
end
tabf = fullfile(outDir, 'segment_table.txt');
fid = fopen(tabf, 'w');  fprintf(fid, '%.13e %.13e %.12e %.9f %.9e %.9e %.12e %.12e %.9e\n', tab');  fclose(fid);

% ---------------- model ----------------
model = ModelUtil.create('T3D');  comp = model.component.create('comp1', true);
try, comp.baseSystem('none'); catch, fb{end+1} = 'baseSystem(none) not accepted: SI units stay on (unit warnings only, same numbers)'; end
em = C.em;  kB = 1.380649e-23;  e0 = 1.602176634e-19;
hconv = 1 / (K.r_convec * L(nl).sx * L(nl).sy);
P = {'D0', em.D0; 'EaJ', em.Ea_eV * e0; 'kBz', kB; 'Bm', em.B; 'Om', em.Omega; 'qe', e0; 'Zs', em.Z_eff; 'rhoCu', em.rho; 'QJ', em.Q_eV * e0; ...
     'Sstep', Sstep; 'Tamb', K.ambient; 'hconv', hconv; 'Tzs', K.T_zero_stress; 'dal', K.dalpha; 'Yr', Y};
for k = 1:size(P, 1), model.param.set(P{k, 1}, sprintf('%.15g', P{k, 2})); end

% ---------------- geometry ----------------
geom = comp.geom.create('geom1', 3);  bn = {};
for i = 1:nl
    nm = sprintf('blk%d', i);  b = geom.create(nm, 'Block');  bn{end+1} = nm; %#ok<AGROW>
    b.set('size', [L(i).sx, L(i).sy, L(i).tm]);
    b.set('pos', [cx - 0.5 * L(i).sx, cy - 0.5 * L(i).sy, zb(i)]);
    if o.inserts
        for q = 1:size(L(i).ins, 1)
            nm = sprintf('ins%d_%d', i, q);  b = geom.create(nm, 'Block');  bn{end+1} = nm; %#ok<AGROW>
            b.set('size', [L(i).ins(q, 2) - L(i).ins(q, 1), L(i).ins(q, 4) - L(i).ins(q, 3), L(i).tm]);
            b.set('pos', [L(i).ins(q, 1), L(i).ins(q, 3), zb(i)]);
        end
    end
end
u1 = geom.create('uni1', 'Union');  u1.selection('input').set(bn);          % one solid object, interior boundaries kept
for q = 1:nr
    a = aj{q};  n = numel(a);
    if isV(q), M = [repmat(x0r(q), n, 1), a, repmat(zr(q), n, 1)]; else, M = [a, repmat(y0r(q), n, 1), repmat(zr(q), n, 1)]; end
    pl = geom.create(sprintf('p%d', q), 'Polygon');  pl.set('type', 'open');
    try
        pl.set('source', 'table');  pl.set('table', M);
    catch
        pl.set('source', 'vectors');  pl.set('x', sprintf('%.15g ', M(:, 1)));  pl.set('y', sprintf('%.15g ', M(:, 2)));  pl.set('z', sprintf('%.15g ', M(:, 3)));
    end
end
geom.feature('fin').set('action', 'assembly');
try, geom.feature('fin').set('createpairs', false); catch, fb{end+1} = 'createpairs not switched off (harmless: no pair is used)'; end
geom.run;
sec_geom = toc(t0);  t0 = tic;

% ---------------- selections by coordinates ----------------
m = 1e-4;  railE = [];  railEq = cell(nr, 1);  dieE = cell(nd, 1);
for d = 1:nd
    i = dieL(d);
    dieE{d} = mphselectbox(model, 'geom1', [cx - 0.5 * L(i).sx - m, cx + 0.5 * L(i).sx + m; cy - 0.5 * L(i).sy - m, cy + 0.5 * L(i).sy + m; ztop(d) + 0.6 * dz, ztop(d) + 2.4 * dz], 'edge');   % strictly between the die top (ztop) and the top of the layer above (>= ztop + 2.5 dz)
    railE = [railE, dieE{d}]; %#ok<AGROW>
end
railE = unique(railE);
if numel(railE) ~= sum(Sr)
    per = arrayfun(@(d) sprintf('%s %d/%d', dieNames{d}, numel(dieE{d}), sum(Sr(die == d - 1))), 1:nd, 'UniformOutput', false);
    error('t3d:sel', 'found %d rail edges, expected %d (one per segment); per die found/expected: %s', numel(railE), sum(Sr), strjoin(per, ', '));
end
ua = sort(unique(across));
if numel(ua) > 1, wsel = 0.2 * min(diff(ua)); else, wsel = 20e-6; end
wsel = max(min(wsel, 20e-6), 1e-6);                           % half-width of the box around one rail
for q = 1:nr
    a = aj{q};
    if isV(q), box = [x0r(q) - wsel, x0r(q) + wsel; a(1) - m, a(end) + m; zr(q) - 0.4 * dz, zr(q) + 0.4 * dz];
    else,      box = [a(1) - m, a(end) + m; y0r(q) - wsel, y0r(q) + wsel; zr(q) - 0.4 * dz, zr(q) + 0.4 * dz]; end
    railEq{q} = mphselectbox(model, 'geom1', box, 'edge');
    if numel(railEq{q}) ~= Sr(q), error('t3d:sel', 'rail %d: found %d edges, expected %d', R(idx(q), 1), numel(railEq{q}), Sr(q)); end
end
dom = cell(1, nl);  topF = cell(1, nl);  botF = cell(1, nl);
for i = 1:nl
    e = 0.05 * min([L.tm]);  bx = [cx - 0.5 * L(i).sx - m, cx + 0.5 * L(i).sx + m; cy - 0.5 * L(i).sy - m, cy + 0.5 * L(i).sy + m];
    dom{i} = mphselectbox(model, 'geom1', [bx; zb(i) - e, zb(i + 1) + e], 'domain');
    botF{i} = mphselectbox(model, 'geom1', [bx; zb(i) - e, zb(i) + e], 'boundary');
    topF{i} = mphselectbox(model, 'geom1', [bx; zb(i + 1) - e, zb(i + 1) + e], 'boundary');
    nexp = 1 + o.inserts * size(L(i).ins, 1);                % at least: the layer + its inserts (an insert that spans the layer splits it)
    if numel(dom{i}) < nexp || (isempty(L(i).ins) && numel(dom{i}) ~= 1) || isempty(botF{i}) || isempty(topF{i})
        error('t3d:sel', 'layer %s: found %d domains (expected %d) / %d bottom faces / %d top faces', L(i).name, numel(dom{i}), nexp, numel(botF{i}), numel(topF{i}));
    end
end
topB = topF{nl};
if numel(topB) ~= 1, error('t3d:sel', 'found %d faces at the top of the sink, expected 1', numel(topB)); end
stackDoms = [];  for i = stackL, stackDoms = [stackDoms, dom{i}(:)']; end %#ok<AGROW>   % dies, bond layers and TIM (what the die-only plots show)

% ---------------- interpolation tables ----------------
names = {'pX0', 'pL', 'pTL', 'pdT', 'pTm', 'plam', 'pj', 'psT'};
try
    f = model.func.create('int1', 'Interpolation');
    f.set('source', 'file');  f.set('filename', tabf);  f.set('nargs', '1');
    for k = 1:numel(names), f.setIndex('funcs', names{k}, k - 1, 0);  f.setIndex('funcs', sprintf('%d', k), k - 1, 1); end
    try, f.set('interp', 'linear');  f.set('extrap', 'const'); catch, end
    try, f.importData; catch, fb{end+1} = 'importData failed: the .mph refers to segment_table.txt instead of embedding it'; end
catch
    try, model.func.remove('int1'); catch, end
    fb{end+1} = 'file-based interpolation not accepted: one table function per quantity';
    xs = arrayfun(@(v) sprintf('%.13e', v), tab(:, 1), 'UniformOutput', false);
    for k = 1:numel(names)
        f = model.func.create(sprintf('int%d', k), 'Interpolation');  f.set('funcname', names{k});
        f.set('table', [xs, arrayfun(@(v) sprintf('%.13e', v), tab(:, k + 1), 'UniformOutput', false)]);
        try, f.set('interp', 'linear');  f.set('extrap', 'const'); catch, end
    end
end
slotExpr = '0';                                               % slot of a rail from (die, direction, coordinate across the rail)
for d = nd:-1:1
    ex = cell(1, 2);
    for kv = 0:1
        g = find(die == d - 1 & isV == kv);  fn = sprintf('pS%d%d', d, kv);
        if isempty(g), t = [0 0; 1 0];
        else, t = [across(g) - 0.7 * W(g), slot(g); across(g) + 0.7 * W(g), slot(g)];  [~, s] = sort(t(:, 1));  t = t(s, :); end
        f = model.func.create(fn, 'Interpolation');  f.set('funcname', fn);
        f.set('table', [arrayfun(@(v) sprintf('%.13e', v), t(:, 1), 'UniformOutput', false), arrayfun(@(v) sprintf('%d', v), t(:, 2), 'UniformOutput', false)]);
        try, f.set('interp', 'linear');  f.set('extrap', 'const'); catch, end
        if kv, ex{2} = [fn '(x)']; else, ex{1} = [fn '(y)']; end
    end
    one = sprintf('if(isV,%s,%s)', ex{2}, ex{1});
    if d == nd, slotExpr = one; else, slotExpr = sprintf('if(z<%.15g,%s,%s)', 0.5 * (ztop(d) + ztop(d + 1)), one, slotExpr); end
end
ztExpr = sprintf('%.15g', ztop(nd));
for d = nd-1:-1:1, ztExpr = sprintf('if(z<%.15g,%.15g,%s)', 0.5 * (ztop(d) + ztop(d + 1)), ztop(d), ztExpr); end

% ---------------- physics 1: stress on the rail lines (dependent variable u, Pa) ----------------
pe = comp.physics.create('ce', 'CoefficientFormEdgePDE', 'geom1');
pe.selection.set(railE);
try, pe.prop('ShapeProperty').set('order', o.order); catch, fb{end+1} = 'element order of the stress not set (COMSOL default)'; end
uname = 'u';  try, uname = char(pe.field('dimensionless').field()); catch, end
% ---------------- HotSpot temperature of every die as a 2-D interpolation function (pictures and tables only) ----------------
Tn = 'Thsv';
hsExpr = '0';
for d = 1:nd
    k = HT(:, 1) == d - 1;  fnm = sprintf('pHT%d', d);
    hf = fullfile(outDir, sprintf('hotspot_T_%s.txt', dieNames{d}));
    fid = fopen(hf, 'w');  fprintf(fid, '%.9e %.9e %.4f\n', HT(k, 2:4)');  fclose(fid);
    try
        f = model.func.create(fnm, 'Interpolation');  f.set('source', 'file');  f.set('filename', hf);  f.set('nargs', '2');
        f.setIndex('funcs', fnm, 0, 0);  f.setIndex('funcs', '1', 0, 1);
        try, f.set('interp', 'linear');  f.set('extrap', 'const'); catch, end
        try, f.importData; catch, end
    catch
        fb{end+1} = sprintf('2-D interpolation of the HotSpot field of %s not created: temperature pictures show 0', dieNames{d});
    end
    if d == nd, hsExpr = sprintf('%s(x,y)', fnm); else, hsExpr = sprintf('if(z<%.15g,%s(x,y),%s)', 0.5 * (ztop(d) + ztop(d + 1)), fnm, hsExpr); end
end
vd = comp.variable.create('var1');  vd.set('Thsv', hsExpr);          % HotSpot temperature painted on the die stack
for fnm = {'pTn'}                                              % filled after the solve (placeholder until then)
    f = model.func.create(fnm{1}, 'Interpolation');  f.set('funcname', fnm{1});  f.set('table', {'0', '300'; '1', '300'});
    try, f.set('interp', 'linear');  f.set('extrap', 'const'); catch, end
end
ve = comp.variable.create('var2');  ve.selection.geom('geom1', 1);  ve.selection.set(railE);     % on the rail lines only
ve.set('ztd',  ztExpr);
ve.set('isV',  sprintf('(z-ztd)>%.15g', 1.5 * dz));
ve.set('alng', 'if(isV,y,x)');
ve.set('Sg',   ['(' slotExpr ')*Sstep+alng']);
ve.set('xloc', 'min(max((alng-pX0(Sg))/pL(Sg),0),1)');
ve.set('axl',  'abs(xloc-0.5)');
ve.set('rch',  'exp((axl-0.5)*plam(Sg))*(1+exp(-2*axl*plam(Sg)))/(1+exp(-plam(Sg)))');
ve.set('rsh',  'sign(xloc-0.5)*exp((axl-0.5)*plam(Sg))*(1-exp(-2*axl*plam(Sg)))/(1+exp(-plam(Sg)))');
ve.set('Tsub', 'pTL(Sg)+pdT(Sg)*xloc');  ve.set('dTsub', 'pdT(Sg)/pL(Sg)');  ve.set('sT0', 'psT(Sg)');
ve.set('Tk',   'Tsub+pTm(Sg)*(1-rch)');
ve.set('dTdx', 'dTsub-pTm(Sg)*plam(Sg)*rsh/pL(Sg)');
ve.set('kap',  'D0*exp(-EaJ/(kBz*Tk))*Bm*Om/(kBz*Tk)');
ve.set('Gd',   'qe*Zs*rhoCu*pj(Sg)/Om');
ve.set('Md',   'QJ/(Om*Tk)*dTdx');
ve.set('tnuc_yr', 'pTn(Sg)');                                  % nucleation time of the rail (years), filled after the solve
ee = pe.feature('cfeq1');
ee.set('c', 'kap');  ee.set('da', 1);  ee.set('ga', {'if(isV,0,kap*(Gd+Md))', 'if(isV,kap*(Gd+Md),0)', '0'});
ee.set('a', 0);  ee.set('f', 0);  ee.set('ea', 0);
pe.feature('init1').set(uname, 'sT0');
sec_build = toc(t0);  t0 = tic;

% ---------------- mesh ----------------
mesh = comp.mesh.create('mesh1');
ed = mesh.create('edg1', 'Edge');  ed.selection.set(railE);                 % rails first: Chebyshev vertices
ds = ed.create('dis1', 'Distribution');  n = o.nper;
try
    ds.set('type', 'explicit');  ds.set('explicit', sprintf('%.12g ', 0.5 * (1 - cos(pi * (0:n) / n))));
catch
    try
        ds.set('type', 'predefined');  ds.set('elemcount', n);  ds.set('elemratio', round(2 * n / pi));  ds.set('method', 'geometric');  ds.set('symmetric', true);
        fb{end+1} = 'explicit distribution not accepted: symmetric geometric grading';
    catch
        ds.set('type', 'number');  ds.set('numelem', n);  fb{end+1} = 'graded distribution not accepted: uniform line elements';
    end
end
meshKind = 'swept prisms';
if truez
    % lateral size per layer, then one sweep per layer from its bottom face upwards (bottom die first, so that the fine
    % triangulation of the die stack is the one that propagates; spreader and sink coarsen outwards)
    for i = 1:nl
        sz = mesh.create(sprintf('size%d', i), 'Size');  sz.selection.geom('geom1', 3);  sz.selection.set(dom{i});
        sz.set('custom', 'on');  sz.set('hmaxactive', true);  sz.set('hmax', sprintf('%.9g', L(i).h));
        sz.set('hminactive', true);  sz.set('hmin', sprintf('%.9g', 0.2 * L(i).h));
    end
    for i = 1:nl
        sw = mesh.create(sprintf('swe%d', i), 'Sweep');  sw.selection.geom('geom1', 3);  sw.selection.set(dom{i});
        try, sw.selection('sourceface').set(botF{i}); catch, fb{end+1} = sprintf('source face of the sweep through %s not set (automatic)', L(i).name); end
        try, sw.selection('targetface').set(topF{i}); catch, end
        dz_ = sw.create('dis1', 'Distribution');  dz_.set('numelem', L(i).nz);
    end
else
    try, mesh.feature('size').set('hauto', o.mesh3d); catch, fb{end+1} = 'predefined mesh size not set (COMSOL default: normal)'; end
    mesh.create('ftet1', 'FreeTet');  meshKind = 'free tetrahedra, scaled bond / TIM thickness';
end
try
    mesh.run;
catch err
    if ~truez, rethrow(err); end
    hint = '';  if ~isempty(IJ), hint = '; or without the inserts: t3d_run(''inserts'',0, ...)'; end
    error('t3d:mesh', ['the swept mesh of the real-thickness stack failed: %s\n' ...
          'Re-run with  t3d_run(''zscale'',''scaled'', ...)  (bond / TIM drawn 50 um thick, free tetrahedra; same temperatures)%s.'], err.message, hint);
end
n_vol = NaN;
try
    ms = mphmeshstats(model, 'mesh1');  n_vol = 0;
    for ty = {'tet', 'prism', 'pyr', 'hex'}
        k_ = find(strcmp(ms.types, ty{1}));  if ~isempty(k_), n_vol = n_vol + double(ms.numelem(k_)); end
    end
catch
    for ty = {'prism', 'tet', 'hex', 'pyr'}
        try, v_ = double(mesh.getNumElem(ty{1}));  if isnan(n_vol), n_vol = 0; end;  n_vol = n_vol + v_; catch, end
    end
end
sec_mesh = toc(t0);  t0 = tic;

sec_thermal = 0;
% ---------------- study 2: transient stress ----------------
std2 = model.study.create('std2');  tt = std2.create('time', 'Transient');
tt.activate('ce', true);
tt.set('tlist', sprintf('%.9g ', tl));  tt.set('rtolactive', true);  tt.set('rtol', sprintf('%g', o.rtol));
try
    sol = model.sol.create('sol2');  sol.study('std2');  sol.attach('std2');  sol.createAutoSequence('std2');
    try, sol.feature('t1').set('initialstepbdfactive', true);  sol.feature('t1').set('initialstepbdf', '1'); catch, fb{end+1} = 'initial BDF step not set (COMSOL default)'; end
    sol.runAll;
catch err
    fb{end+1} = ['solver sequence of the stress study not customised: ' err.message];
    try, model.sol.remove('sol2'); catch, end
    std2.run;
end
sec_stress = toc(t0);  t0 = tic;
tg = cell(model.result.dataset.tags());  dsU = tg{end};
for k = numel(tg):-1:1
    try, info = mphsolinfo(model, 'dataset', tg{k});  if numel(info.solvals) == nt, dsU = tg{k};  break; end, catch, end
end
dsT = dsU;

% ---------------- results: temperature ----------------
Tc = HT(:, 4);                                               % no thermal solve: both columns are the HotSpot field
fid = fopen(fullfile(outDir, 't3d_die_temperature.csv'), 'w');  fprintf(fid, 'die_id,x_m,y_m,T_comsol_K,T_hotspot_K\n');
fprintf(fid, '%d,%.9e,%.9e,%.4f,%.4f\n', [HT(:, 1:3), Tc(:), HT(:, 4)]');  fclose(fid);
for d = 1:nd
    k = HT(:, 1) == d - 1;
    fprintf('%s temperature (HotSpot, used as is): %.2f .. %.2f K\n', dieNames{d}, min(HT(k, 4)), max(HT(k, 4)));
end
% ---------------- results: the whole line field at every output time (time x position) ----------------
no = [0; cumsum(Sr + 1)];  X = zeros(no(end), 3);  Ths = zeros(no(end), 1);
for q = 1:nr
    k = no(q)+1:no(q+1);  a = aj{q};  rows = R(idx(q), 14) + (1:Sr(q));
    if isV(q), X(k, :) = [repmat(x0r(q), numel(a), 1), a, repmat(zr(q), numel(a), 1)]; else, X(k, :) = [a, repmat(y0r(q), numel(a), 1), repmat(zr(q), numel(a), 1)]; end
    Ths(k) = [S(rows, 5); S(rows(end), 6)];
end
Tnode = Ths;
PTS = [];  FLD = [];  loc = [];
for i0 = 1:10:nt
    sn = i0:min(nt, i0 + 9);
    dd = mpheval(model, uname, 'edim', 1, 'selection', railE, 'dataset', dsU, 'solnum', sn);
    if isempty(PTS)
        PTS = dd.p;  [~, u_] = unique(round(PTS' * 1e9), 'rows', 'stable');  PTS = PTS(:, u_);  FLD = zeros(nt, size(PTS, 2), 'single');
        [ok, loc] = ismember(round(X * 1e9), round(PTS' * 1e9), 'rows');
        if ~all(ok), error('t3d:nodes', '%d rail nodes were not found among the mesh vertices', nnz(~ok)); end
    end
    v = dd.d1;  if size(v, 1) ~= numel(sn), v = v.'; end
    FLD(sn, :) = single(v(:, u_));
end
sig = double(FLD(:, loc));                                    % nt x nodes
% rail id and position along the rail of every field point
prail = zeros(1, size(PTS, 2));  ppos = zeros(1, size(PTS, 2));
for q = 1:nr
    if isV(q), mq = abs(PTS(3, :) - zr(q)) < 0.25 * dz & abs(PTS(1, :) - x0r(q)) < 1e-9 & PTS(2, :) >= aj{q}(1) - 1e-9 & PTS(2, :) <= aj{q}(end) + 1e-9;
    else,      mq = abs(PTS(3, :) - zr(q)) < 0.25 * dz & abs(PTS(2, :) - y0r(q)) < 1e-9 & PTS(1, :) >= aj{q}(1) - 1e-9 & PTS(1, :) <= aj{q}(end) + 1e-9; end
    prail(mq) = R(idx(q), 1);  if isV(q), ppos(mq) = PTS(2, mq) - aj{q}(1); else, ppos(mq) = PTS(1, mq) - aj{q}(1); end
end
% ---------------- consistency checks ----------------
if useHS, Tchk = Ths; else, Tchk = Tnode; end
u0 = max(C.em.B * K.dalpha * (K.T_zero_stress - Tchk'), 0);
if useHS, for q = 1:nr, rows = R(idx(q), 14) + (1:Sr(q));  k = no(q)+1:no(q+1);  u0(k) = NaN;  u0(k(2:end-1)) = 0.5 * (S(rows(1:end-1), 9) + S(rows(2:end), 9))'; end, end
du0 = max(abs(sig(1, :) - u0), [], 'omitnan');
fprintf('initial stress at the rail nodes: %.1f .. %.1f MPa (largest deviation from the expected thermal residual stress %.2f MPa)\n', min(sig(1, :)) / 1e6, max(sig(1, :)) / 1e6, du0 / 1e6);
if du0 > 5e6, error('t3d:init', 'the initial stress differs from the thermal residual stress by %.1f MPa: the rails did not get the temperature of the thermal solution', du0 / 1e6); end
try
    % checked at the segment mid points (a junction sits in the 4-nm ramp between two table rows, where the per-segment
    % table functions are not meant to be read; no integration point ever lies there)
    XM = zeros(3, sum(Sr));  Texp = zeros(1, sum(Sr));  p = 0;
    for q = 1:nr
        a = aj{q};  rows = R(idx(q), 14) + (1:Sr(q));  am = 0.5 * (a(1:end-1) + a(2:end))';  k = p + (1:Sr(q));  p = p + Sr(q);
        if isV(q), XM(:, k) = [repmat(x0r(q), 1, Sr(q)); am; repmat(zr(q), 1, Sr(q))]; else, XM(:, k) = [am; repmat(y0r(q), 1, Sr(q)); repmat(zr(q), 1, Sr(q))]; end
        Texp(k) = 0.5 * (S(rows, 5) + S(rows, 6))';
    end
    v = mphinterp(model, 'Tsub', 'coord', XM, 'dataset', dsU, 'solnum', nt, 'edim', 1);  v = v(:)';  dTk = max(abs(v - Texp));
    fprintf('rail temperature inside the stress model against its source at the segment mid points: largest difference %.3f K\n', dTk);
    if dTk > 0.2, error('t3d:Tk', 'the stress model uses a rail temperature that differs from its source by %.2f K', dTk); end
catch err
    if strcmp(err.identifier, 't3d:Tk'), rethrow(err); end
    fb{end+1} = 'rail temperature inside the stress model could not be read back for the check';
end
% ---------------- results: nucleation, per-rail and per-node tables ----------------
t_nuc = inf(nr, 1);
fr = fopen(fullfile(outDir, 't3d_rail_result.csv'), 'w');  fprintf(fr, 'rail,die_id,kind_is_V,index,n_seg,t_nuc_3d_s');
for h = 1:nh, fprintf(fr, ',smax_node_%gyr_Pa', hz(h)); end;  fprintf(fr, '\n');
fs = fopen(fullfile(outDir, 't3d_rail_smax.csv'), 'w');  fprintf(fs, 'rail');  fprintf(fs, ',%.9g', tl);  fprintf(fs, '\n');
fn = fopen(fullfile(outDir, 't3d_nodes.csv'), 'w');  fprintf(fn, 'rail,node,die_id,x_m,y_m,T_3d_K,T_hotspot_K');
for h = 1:nh, fprintf(fn, ',sigma_%gyr_Pa', hz(h)); end;  fprintf(fn, '\n');
nid = zeros(1, no(end));
for q = 1:nr
    k = no(q)+1:no(q+1);  sg = sig(:, k);  smax = max(sg, [], 2);  kc = find(smax >= sc, 1);  nid(k) = R(idx(q), 1) * 1000 + (0:Sr(q));
    if ~isempty(kc)
        if kc == 1, t_nuc(q) = tl(1);
        elseif smax(kc) > smax(kc-1), t_nuc(q) = (sqrt(tl(kc-1)) + (sc - smax(kc-1)) / (smax(kc) - smax(kc-1)) * (sqrt(tl(kc)) - sqrt(tl(kc-1))))^2;
        else, t_nuc(q) = tl(kc); end
    end
    fprintf(fr, '%d,%d,%d,%d,%d,%.9e', R(idx(q), 1), die(q), isV(q), R(idx(q), 4), Sr(q), t_nuc(q));  fprintf(fr, ',%.6e', smax(ih));  fprintf(fr, '\n');
    fprintf(fs, '%d', R(idx(q), 1));  fprintf(fs, ',%.6e', smax);  fprintf(fs, '\n');
    fprintf(fn, ['%d,%d,%d,%.9e,%.9e,%.4f,%.4f' repmat(',%.6e', 1, nh) '\n'], ...
            [repmat(R(idx(q), 1), 1, numel(k)); 0:Sr(q); repmat(die(q), 1, numel(k)); X(k, 1)'; X(k, 2)'; Tnode(k)'; Ths(k)'; sg(ih, :)]);
end
fclose(fr);  fclose(fs);  fclose(fn);
% time x node history (csv) and the full time x position field (mat)
fid = fopen(fullfile(outDir, 't3d_nodes_history.csv'), 'w');
fprintf(fid, 'time_s');  fprintf(fid, ',%d', nid);  fprintf(fid, '\n');          % column id = rail*1000 + node
for it = 1:nt, fprintf(fid, '%.9g', tl(it));  fprintf(fid, ',%.6e', sig(it, :));  fprintf(fid, '\n'); end
fclose(fid);
if o.history
    hist = struct('case_name', caseName, 'times_s', tl, 'horizons_years', hz, 'node_id', nid, 'node_xyz', X, 'node_rail', floor(nid / 1000), 'node_T3d_K', Tnode', 'node_Thotspot_K', Ths', ...
                  'sigma_nodes_Pa', single(sig), 'pt_xyz', PTS', 'pt_rail', prail, 'pt_pos_m', ppos, 'sigma_points_Pa', FLD, ...
                  'rail_id', R(idx, 1)', 't_nuc_s', t_nuc', 'note', 'sigma_points_Pa(it, ip): stress at output time it at line-mesh vertex ip (60 per segment); sigma_nodes_Pa: at the segment end points (junctions)'); %#ok<NASGU>
    save(fullfile(outDir, 't3d_history.mat'), '-struct', 'hist', '-v7.3');
end
% field along every segment at the horizons (same format as the 1-D / 2-D runs)
for h = 1:nh
    fid = fopen(fullfile(outDir, sprintf('t3d_field_%gyr.csv', hz(h))), 'w');  fprintf(fid, 'rail,seg,xloc,s_m,sigma_Pa\n');
    xi = linspace(0, 1, o.field_pts);
    for q = 1:nr
        a = aj{q};  mq = prail == R(idx(q), 1);  pp = ppos(mq);  vv = double(FLD(ih(h), mq));  [pp, s_] = sort(pp);  vv = vv(s_);
        for k = 1:Sr(q)
            sloc = (a(k) - a(1)) + xi * (a(k + 1) - a(k));  val = interp1(pp, vv, sloc, 'linear', 'extrap');
            fprintf(fid, '%d,%d,%.4f,%.9e,%.6e\n', [repmat(R(idx(q), 1), 1, o.field_pts); repmat(k - 1, 1, o.field_pts); xi; sloc; val]);
        end
    end
    fclose(fid);
end
% nucleation time as a function on the rails (for the COMSOL plot)
tn_tab = zeros(2 * nr, 2);
for q = 1:nr
    a = aj{q};  v = t_nuc(q) / Y;  if ~isfinite(v), v = 10 * o.tend_years; end
    tn_tab(2 * q - 1, :) = [slot(q) * Sstep + a(1) + ep, v];  tn_tab(2 * q, :) = [slot(q) * Sstep + a(end) - ep, v];
end
try
    model.func('pTn').set('table', [arrayfun(@(v) sprintf('%.13e', v), tn_tab(:, 1), 'UniformOutput', false), arrayfun(@(v) sprintf('%.6g', v), tn_tab(:, 2), 'UniformOutput', false)]);
catch
    fb{end+1} = 'nucleation-time table not filled (the COMSOL nucleation plot shows the placeholder)';
end
sec_post = toc(t0);  t0 = tic;

% ---------------- plot groups and COMSOL-rendered images ----------------
pgs = {};
solT = '';  solU = '';
try, solU = char(model.result.dataset(dsU).getString('solution'));  solT = solU; catch, end
% data sets restricted to the die stack (temperature) and to the rails (stress): what the plots show, and what the
% automatic view frames.  The spreader and sink stay in the model; they only leave the pictures.
dsTd = dsT;  dsUr = dsU;
try
    d_ = model.result.dataset.create('dsetStack', 'Solution');  d_.set('solution', solT);
    d_.selection.geom('geom1', 3);  d_.selection.set(stackDoms);  dsTd = 'dsetStack';
    d_ = model.result.dataset.create('dsetRails', 'Solution');  d_.set('solution', solU);
    d_.selection.geom('geom1', 1);  d_.selection.set(railE);  dsUr = 'dsetRails';
    for d = 1:nd
        d_ = model.result.dataset.create(sprintf('dsetRails%d', d), 'Solution');  d_.set('solution', solU);
        d_.selection.geom('geom1', 1);  d_.selection.set(dieE{d});
    end
catch
    fb{end+1} = 'selection data sets not created: the COMSOL pictures show the whole package (zoom in the .mph)';
end
vtop = '';                                                    % an orthographic view from above for the per-die pictures
if o.topview
    try
        vw = model.view.create('vtop', 3);  Sd = max(SX, SY);  Dcam = 10 * Sd;  zt_ = zb(dieL(nd) + 1);
        vw.camera.set('projection', 'orthographic');
        vw.camera.set('position', [cx, cy, zt_ + Dcam]);  vw.camera.set('target', [cx, cy, zt_]);  vw.camera.set('up', [0 1 0]);
        vw.camera.set('rotationpoint', [cx, cy, zt_]);  vw.camera.set('viewoffset', [0 0]);
        vw.camera.set('zoomangle', 2 * atand(0.62 * Sd / Dcam));          % frames the die with a small margin
        vtop = 'vtop';
    catch
        fb{end+1} = 'top view not created: per-die pictures use the default 3-D view';
    end
end
try   % temperature of the whole package
    pg = model.result.create('pgT_all', 'PlotGroup3D');  pg.label('HotSpot temperature on the package (K; package faces show 0)');  pg.set('data', dsT);
    sf = pg.create('surf1', 'Surface');  sf.set('expr', Tn);  sf.set('colortable', 'Inferno');  pgs{end+1} = 'pgT_all';
catch, fb{end+1} = 'plot: temperature of the package not created'; end
try   % temperature of the die stack only
    pg = model.result.create('pgT_stack', 'PlotGroup3D');  pg.label('HotSpot temperature painted on the die stack (K)');  pg.set('data', dsTd);
    sf = pg.create('surf1', 'Surface');  sf.set('expr', Tn);  sf.set('colortable', 'Inferno');  pgs{end+1} = 'pgT_stack';
catch, fb{end+1} = 'plot: temperature of the die stack not created'; end
tube3 = sprintf('%.3g', 1.5e-3 * max(SX, SY));  tube2 = sprintf('%.3g', 2e-3 * max(SX, SY));     % tube radius of the rail plots (6 / 8 um on a 4-mm die)
for h = 1:nh   % 3-D stress on the rails above grey dies, one picture per horizon
    try
        tgp = sprintf('pgS3d_%d', h);  pg = model.result.create(tgp, 'PlotGroup3D');  pg.label(sprintf('Rail stress after %g years, 3-D (Pa)', hz(h)));
        pg.set('data', dsUr);  pg.setIndex('looplevel', ih(h), 0);
        ln = pg.create('line1', 'Line');  ln.set('expr', uname);  ln.set('colortable', 'RainbowLight');
        try, ln.set('linetype', 'tube');  ln.set('radiusexpr', tube3);  ln.set('tuberadiusscaleactive', true);  ln.set('tuberadiusscale', 1); catch, end
        try, ln.set('rangecoloractive', true);  ln.set('rangecolormin', -sc);  ln.set('rangecolormax', sc); catch, end
        try
            sd = pg.create('surf1', 'Surface');  sd.set('data', dsTd);  sd.set('expr', '1');  sd.set('coloring', 'uniform');  sd.set('color', 'gray');
            sd.set('alpha', 0.25);  try, sd.set('solutionparams', 'parent'); catch, end
        catch
        end
        pgs{end+1} = tgp;
    catch, fb{end+1} = sprintf('plot: 3-D rail stress after %g years not created', hz(h)); end
end
for d = 1:nd   % per-die top view at every horizon (the 2-D stress map of that die's power grid)
    for h = 1:nh
        try
            tgp = sprintf('pgS2d_%s_%d', dieNames{d}, h);  pg = model.result.create(tgp, 'PlotGroup3D');
            pg.label(sprintf('%s rail stress after %g years, top view (Pa)', dieNames{d}, hz(h)));
            pg.set('data', sprintf('dsetRails%d', d));  pg.setIndex('looplevel', ih(h), 0);
            if ~isempty(vtop), pg.set('view', vtop); end
            ln = pg.create('line1', 'Line');  ln.set('expr', uname);  ln.set('colortable', 'RainbowLight');
            try, ln.set('linetype', 'tube');  ln.set('radiusexpr', tube2);  ln.set('tuberadiusscaleactive', true);  ln.set('tuberadiusscale', 1); catch, end
            try, ln.set('rangecoloractive', true);  ln.set('rangecolormin', -sc);  ln.set('rangecolormax', sc); catch, end
            pgs{end+1} = tgp;
        catch, fb{end+1} = sprintf('plot: top view of %s after %g years not created', dieNames{d}, hz(h)); end
    end
end
try   % nucleation time of every rail, 3-D
    pg = model.result.create('pgTn', 'PlotGroup3D');  pg.label('Nucleation time of every rail (years, log10)');  pg.set('data', dsUr);  pg.setIndex('looplevel', nt, 0);
    ln = pg.create('line1', 'Line');  ln.set('expr', 'log10(tnuc_yr)');  ln.set('colortable', 'Viridis');  ln.set('colortablerev', true);
    try, ln.set('linetype', 'tube');  ln.set('radiusexpr', tube3);  ln.set('tuberadiusscaleactive', true);  ln.set('tuberadiusscale', 1); catch, end
    try, ln.set('rangecoloractive', true);  ln.set('rangecolormin', log10(0.1));  ln.set('rangecolormax', log10(10 * o.tend_years)); catch, end
    pgs{end+1} = 'pgTn';
catch, fb{end+1} = 'plot: nucleation time not created'; end
nimg = 0;
for k = 1:numel(pgs)
    try
        im = model.result.export.create(sprintf('img%d', k), 'Image');
        im.set('sourceobject', pgs{k});  im.set('filename', fullfile(imgDir, [pgs{k} '.png']));
        try, im.set('size', 'manualweb');  im.set('unit', 'px');  im.set('width', '2400');  im.set('height', '1800'); catch, end
        im.run;  nimg = nimg + 1;
    catch
        fb{end+1} = ['COMSOL image of plot group ' pgs{k} ' not exported (open the .mph and export it from the plot window)'];
    end
end
if o.rail_plots   % one 1-D stress plot per rail, all horizons as curves
    railImg = fullfile(imgDir, 'rails');  if ~exist(railImg, 'dir'), mkdir(railImg); end
    nfail = 0;
    for q = 1:nr
        try
            tgp = sprintf('pgR%d', q);  pg = model.result.create(tgp, 'PlotGroup1D');
            pg.label(sprintf('rail %d (%s %s%d): stress along the rail', R(idx(q), 1), dieNames{die(q) + 1}, char('H' + isV(q) * ('V' - 'H')), R(idx(q), 4)));
            pg.set('data', dsU);  pg.set('looplevelinput', {'manual'});  pg.set('looplevel', {sprintf('%d ', ih)});
            lg = pg.create('lngr1', 'LineGraph');  lg.selection.set(railEq{q});  lg.set('expr', uname);
            lg.set('xdata', 'expr');  if isV(q), lg.set('xdataexpr', 'y'); else, lg.set('xdataexpr', 'x'); end
            try, lg.set('legend', true);  lg.set('legendmethod', 'automatic'); catch, end
            try, pg.set('xlabel', 'position along the rail (m)');  pg.set('ylabel', 'hydrostatic stress (Pa)');  pg.set('xlabelactive', true);  pg.set('ylabelactive', true); catch, end
            if q == 1, pgs{end+1} = tgp; end
            im = model.result.export.create(sprintf('imgR%d', q), 'Image');
            im.set('sourceobject', tgp);  im.set('filename', fullfile(railImg, sprintf('rail_%05d_%s_%s%d.png', R(idx(q), 1), dieNames{die(q) + 1}, char('H' + isV(q) * ('V' - 'H')), R(idx(q), 4))));
            try, im.set('size', 'manualweb');  im.set('unit', 'px');  im.set('width', '1600');  im.set('height', '900'); catch, end
            im.run;  nimg = nimg + 1;
            if q > 1, model.result.remove(tgp);  model.result.export.remove(sprintf('imgR%d', q)); end      % keep one in the .mph as a template
        catch
            nfail = nfail + 1;
        end
    end
    if nfail, fb{end+1} = sprintf('%d of %d per-rail plots could not be exported', nfail, nr); end
end
try, model.label(sprintf('%s stack: rail stress with HotSpot temperatures (no thermal solve)', caseName)); catch, end
mphsave(model, fullfile(outDir, sprintf('%s_3d.mph', caseName)));
sec_plots = toc(t0);

meta = struct('case_name', caseName, 'options', o, 'n_rails', nr, 'n_segments', sum(Sr), 'n_line_elements', sum(Sr) * o.nper, 'n_volume_elements', n_vol, 'volume_mesh', meshKind, ...
              'times_s', tl, 'horizons_years', hz, 'n_images', nimg, ...
              'seconds', struct('geometry', sec_geom, 'build', sec_build, 'mesh', sec_mesh, 'thermal', sec_thermal, 'stress', sec_stress, 'post', sec_post, 'plots_and_save', sec_plots), ...
              'temperature_source', 'HotSpot (no thermal solve in COMSOL)', 'stress_variable', uname, 'z_layers_model_m', zb, 'layer_names', {{L.name}}, 'layer_kinds', {{L.kind}}, 'n_inserts', numel(IJ), 'fallbacks', {fb});
try, meta.comsol = char(ModelUtil.getComsolVersion()); catch, end
fid = fopen(fullfile(outDir, 't3d_meta.json'), 'w');  fprintf(fid, '%s', jsonencode(meta));  fclose(fid);
fprintf('3-D model of %s (HotSpot temperatures): %d rails, %d line elements, %g volume elements for the pictures (%s) | geometry %.1f build %.1f mesh %.1f | STRESS %.1f s (%d output times to %g years) | post %.1f plots %.1f\n', ...
        caseName, nr, sum(Sr) * o.nper, n_vol, meshKind, sec_geom, sec_build, sec_mesh, sec_stress, nt, o.tend_years, sec_post, sec_plots);
fprintf('earliest nucleation %.3f years; %d rails nucleate within %g years; %d pictures in %s\n', min(t_nuc) / Y, nnz(t_nuc <= o.tend_years * Y), o.tend_years, nimg, imgDir);
if ~isempty(fb), fprintf('API fall-backs taken:\n');  fprintf('   %s\n', fb{:}); end
fprintf('results in %s\n', outDir);
out = meta;
end
