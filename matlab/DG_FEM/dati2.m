%=======================================================================================================
% This contain all the information for running main
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================
%
%  DATI= struct( 'name',              % set the name of the test  
%                'method',            % (string) e.g. 'SIP','NIP'or 'IIP'
%                'Domain',            % set the domain [x1,x2;y1,y2]
%                'exact_sol',         % set the exact solution
%                'source',            % set the forcing term
%                'grad_exact_1',      % set the first componenet of the gradient of the exact solution
%                'grad_exact_2',      % set the second componenet of the gradient of the exact solution
%                'fem',               % set finite element space (e.g.,
%                'P1', 'P2', 'P3')
%                'penalty_coeff'      % (real) penalty pearameter
%                'nqn',               % (integer) number of 1D Gauss-Ledendre quadrature nodes in 1 
%                       dimension [1,2,3,4,5,6,.....]
%                'nqn_2D',          % number of quadrature nodes for integrals over elements
%========================================================================================================

function [DATA] = dati(test)

if test=='Test1'
    
DATA = struct(  'name',             test,...
               ... % Test name
                'method',           'NIP',...  
               ... % Set DG discretization
               'domain',           [0,1;0,1],...
               ... % Domain bounds
               'exact_sol',        '(x-x.^2).*exp(3*x).*sin(2*pi*y)',...
               ... % Definition of exact solution
               'source',           '(-4 + (3+4*pi^2)*x + (9-4*pi^2)*x.^2).*exp(3*x).*sin(2*pi*y)',...
               ... % Forcing term
               'grad_exact_1',     '(1+x-3*x.^2).*exp(3*x).*sin(2*pi*y)',... 
               ... % Definition of exact gradient (x comp) 
               'grad_exact_2',     '2*pi*(x-x.^2).*exp(3*x).*cos(2*pi*y)',...    
               ... % Definition of exact gradient (y comp)
               'fem',              'P2',...   
               ... % Finite element space (other choices 'P2', 'P3')'
               'penalty_coeff',     1,... 
               ... % Penalty coefficient
               'nqn',               4 ...
               ... % Number of 1d GL quadrature nodes
               );
           
end
