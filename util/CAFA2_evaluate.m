addpath(genpath('..\analysis\CAFA2\github\CAFA2-master'))
cd('E:\swang141\project\SequencingNetwork\Sheng\src');
config = cafa_parse_config('example.job');
cafa_driver_preeval(config)
cafa_driver_eval(config)