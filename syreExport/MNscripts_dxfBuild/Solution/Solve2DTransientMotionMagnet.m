function magDataUpd = Solve2DTransientMotionMagnet(magData)

% Solve2DStaticMagnet.m [v1.00.00 (30-11-2012)]
% 2D magnetostatic solution
% =========================================================================
% Syntax: magDataUpd = Solve2DStaticMagnet(magData)
% Input:
%         - magData:    Magnet's data structure
% Output:
%         - magDataUpd: updated Magnet's data structure
% =========================================================================

dh = magData.documentHandler;
invoke(dh,'solveTransient2dWithMotion');
sh = invoke(dh, 'getSolution');
magDataUpd = magData;
magDataUpd.solutionHandler = sh;