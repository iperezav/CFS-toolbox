function varargout=CFS(varargin)

%CFS Chen-Fliess series of a nonlinear affine control system 
%   varargout{1} = CFS(h,g,z,N) provides the collection of coefficients from 
%   the Chen-Fliess series, excluding the coefficient (c, \emptyset), 
%   for the nonlinear affine control system specified by
%               \dot{z}=g_0(z)+sum_{i=1}^m g_i(z) u_i
%                     y=h(z)
%                  z(0)=z_0
%   with g=[g_0 ... g_m] and N is the truncation length of the words of 
%   the Chen-Fliess series.
%
%   varargout{1} = CFS(u,T,N) returns the list of iterative integrals up to 
%   words of length N
%    E_{x_i\eta}[u](t_0,t)=int_{t0}^t u_i(\tau)E_{\eta}[u](t_0,\tau)d\tau_2
%   for all |eta|<=N and t<=T where u=[u_0 u_1 ... u_m], u_0= 1 is added in
%   the code and only u_1,...,u_m are input arguments of the function.
%
%   [varargout{1},varargout{2}] = CFS(h,g,z,N,u,t) returns both the list
%   coefficients of the nonlinear affine control system and the list of 
%   iterated integrals of words up to length N.
%
%   Example 1:
%      syms x [1 2];
%      g0 = [-x1*x2;x1*x2];
%      g1 = [x1;0];
%      g2 = [0;-x2];
%      h = x1;
%      g = [g0 g1 g2];
%      N = 2;
%      coeff = CFS(h,g,z,N);
%      returns the 12-by-1 symbolic array
%      where each row is (c,\eta) for \eta from x0 to x2x2
%                      -x1*x2
%                          x1
%                           0
%         - x1^2*x2 + x1*x2^2
%                      -x1*x2
%                       x1*x2
%                      -x1*x2
%                          x1
%                           0
%                           0
%                           0
%                           0
%
%   Example 2:
%      dt = 0.001;
%      t0 = 0;
%      tf =3 ;
%      t = t0:dt:tf;
%      N = 2;
%      u1 = sin(t);
%      u2 = cos(t);
%      u = [u1;u2];
%      intt = CFS(u,t,N);
%      returns a 12-by-3001 array 
%      where each row is E_{\eta}[u](t) for a fixed \eta from x0 to x2x2 for all
%      t \in [t0 tf] and each column represents E_{\eta}[u](T) for 
%      a fixed T \in [t0 tf] and all \etas from x0 to x2x2. Here, only a
%      sample of 12-by-5 matrix is shown
%         0	-0.00100000000000000	-0.00200000000000000	-0.00300000000000000	-0.00400000000000000 
%         0	-0.00100000000000000	-0.00200000000000000	-0.00300000000000000	-0.00400000000000000
%         0	-0.00100000000000000	-0.00200000000000000	-0.00300000000000000	-0.00400000000000000
%         0	0	1.00000000000000e-06	3.00000000000000e-06	6.00000000000000e-06
%         0	0	1.00000000000000e-06	3.00000000000000e-06	6.00000000000000e-06
%         0	0	1.00000000000000e-06	3.00000000000000e-06	6.00000000000000e-06
%         0	0	1.00000000000000e-06	3.00000000000000e-06	6.00000000000000e-06
%         0	0	1.00000000000000e-06	3.00000000000000e-06	6.00000000000000e-06
%         0	0	1.00000000000000e-06	3.00000000000000e-06	6.00000000000000e-06
%         0	0	1.00000000000000e-06	3.00000000000000e-06	6.00000000000000e-06
%         0	0	1.00000000000000e-06	3.00000000000000e-06	6.00000000000000e-06
%         0	0	1.00000000000000e-06	3.00000000000000e-06	6.00000000000000e-06
%
%   Example 3:
%      syms x [1 2];
%      g0 = [-x1*x2;x1*x2];
%      g1 = [x1;0];
%      g2 = [0;-x2];
%      h = x1;
%      g = [g0 g1 g2];
%      dt = 0.001;
%      t0 = 0;
%      tf = 3;
%      t = t0:dt:tf;
%      N = 2;
%      u1 = sin(t);
%      u2 = cos(t);
%      u = [u1;u2];
%      [coeff intt]=CFS(h,g,z,N,u,t);
%      returns the previous results in Example 1 and 2 in coeff and intt
%      respectively.


narginchk(3,6);

if nargin==4
    h_output=varargin{1};
    v_field=varargin{2};
    v_states=varargin{3};
    N_trunc=varargin{4};
    sch_c=sanity_check_coef(h_output,v_field,v_states,N_trunc);
    if sch_c

        varargout{1}=power_series(h_output,v_field,v_states,N_trunc);
    else
        error(message('Error computing the coefficients'));
    end
elseif nargin==3
    v_u=varargin{1};
    v_t=varargin{2};
    N_trunc=varargin{3};
    sch_i=sanity_check_itint(v_u,v_t,N_trunc);
    if sch_i
        varargout{1} = iterative_int(v_u,v_t,N_trunc);
    else
        error(message('Error computing the iterative integrals'));
    end
elseif nargin==6
    h_output=varargin{1};
    v_field=varargin{2};
    v_states=varargin{3};
    N_trunc=varargin{4};
    v_u=varargin{5};
    v_t=varargin{6};
    sch_ci=sanity_check_coefit(h_output,v_field,v_states,N_trunc,v_u,v_t);
    if sch_ci
        varargout{1}=power_series(h_output,v_field,v_states,N_trunc);
        varargout{2} = iterative_int(v_u,v_t,N_trunc);
    else 
        error(message('Error computing the iterative integrals and the coefficients'));
    end
elseif nargin==5
    error(message('Too many arguments'));

end


    function sch_c=sanity_check_coef(h_output,v_field,v_states,N_trunc)
    
        sch_c=isa(h_output, 'sym') & isa(v_field, 'sym') & isa(v_states, 'sym') &...
                size(v_states,2)==size(v_field,1) & N_trunc>0 &...
                floor(N_trunc)== N_trunc ;
    
    end
    
    
    function sch_i=sanity_check_itint(v_u,v_t,N_trunc)
    
        sch_i=~isa(v_t, 'sym') & ~isa(N_trunc, 'sym') & ~isa(v_u, 'sym') &...
            size(v_u,2)==size(v_t,2) & N_trunc>0 & floor(N_trunc)== N_trunc &...
                size(v_t,1)==1;    
    end

    function sch_ci=sanity_check_coefit(h_output,v_field,v_states,N_trunc,v_u,v_t)
    
        sch_ci=sanity_check_coef(h_output,v_field,v_states,N_trunc);
        sch_ci=sch_ci*sanity_check_itint(v_u,v_t,N_trunc);


    
    end



end