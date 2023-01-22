function conf = CRC_init(conf)

if ~isfield(conf,'MaxIter'), conf.MaxIter = 500; end;
if ~isfield(conf,'gamma'), conf.gamma = 0.9; end;
if ~isfield(conf,'lambda'), conf.lambda = 1.0; end; % 1.0
if ~isfield(conf,'theta'), conf.theta = 0.75; end;
if ~isfield(conf,'a'), conf.a = 1.0; end;
if ~isfield(conf,'ecr'), conf.ecr = 1e-5; end;
if ~isfield(conf,'minP'), conf.minP = 1e-5; end;
