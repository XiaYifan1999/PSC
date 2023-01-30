function conf = CRC_LLT_init(conf)

if ~isfield(conf,'MaxIter'), conf.MaxIter = 500; end;
if ~isfield(conf,'gamma'), conf.gamma = 0.9; end;
if ~isfield(conf,'lambda1'), conf.lambda1 = 1; end;
if ~isfield(conf,'lambda2'), conf.lambda2 = 1000;end; 
if ~isfield(conf,'tau'), conf.tau = 0.75; end;
if ~isfield(conf,'a'), conf.a = 1.0; end;
if ~isfield(conf,'ecr'), conf.ecr = 1e-5; end;
if ~isfield(conf,'minP'), conf.minP = 1e-5; end;
if ~isfield(conf,'Kn'), conf.Kn = 15; end;