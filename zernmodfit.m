function [ad,varargout] = zernmodfit(r,theta,data,N)
%ZERNMODFIT Fit data to a modified Zernike basis.
%   AD = ZERNMODFIT(R,THETA,DATA,N) returns a two-column matrix containing 
%   the modal coefficients A and rotation angles (orientation axes) D for
%   DATA expressed at the coordinate locations (R,THETA) within the unit 
%   disk.  The A and D are computed for all [n m] modes with n = 0 to N, 
%   and m = 0:2:n for even n, m = 1:2:n for odd n.  R, THETA, and DATA must
%   be vectors, and must contain the same number of elements.  N must be a 
%   positive integer.
% 
%   [AD,NM] = ZERNMODFIT(R,THETA,DATA,N) returns a 2-column matrix NM
%   containing the [n m] mode numbers associated with each (coefficient,
%   angle) pair in the array AD.
% 
%   Note: For the m > 0 modes, the computed Zernike coefficients A are
%   always positive, since the orientation axis D for each mode can always
%   be chosen so as to make them positive.  For the m = 0 modes, which
%   have no angular dependence, the coefficients can be both positive and 
%   negative.
% 
%   Example
% 
%     % Create some data
%     N = 100;
%     x = (-N:2:N)/N;
%     z = peaks(N+1);
%     [X,Y] = meshgrid(x);
%     [theta,r] = cart2pol(X,Y);
%     is_in = r <= 1;
%     z(~is_in) = nan;
%     r = r(is_in);
%     theta = theta(is_in);
% 
%     % Compute the modified Zernike fit
%     N = 9;
%     [ad,nm] = zernmodfit(r,theta,z(is_in),N);
% 
%     % Reconstruct the data
%     zr = z; zr(is_in) = 0;
%     Nz = size(ad,1);
%     for k = 1:Nz
%         zr(is_in) = ...
%         zr(is_in) + ad(k,1)*zernfun(nm(k,1),nm(k,2),r,theta + ad(k,2));
%     end
% 
%     % Display the original and reconstructed functions
%     h = figure('Visible','off','Position',[0 0 1000 500]);
%     subplot(1,2,1)
%     pcolor(x,x,z), shading interp
%     set(gca,'XTick',[],'YTick',[])
%     axis square
%     title('Original')
% 
%     subplot(1,2,2)
%     pcolor(x,x,zr), shading interp
%     set(gca,'XTick',[],'YTick',[])
%     axis square
%     title(['Reconstructed (N = ' num2str(N) ')'])
%     movegui(h,'center'), set(h,'Visible','on')
% 
%     % Display the modified Zernike spectrum
%     h = figure('Visible','off','Position',[0 0 1200 400]);
%     stem(ad(:,1))
%     set(gca,'TickLength',[0 0],'XTick',[],'XLim',[0 Nz+1],'YGrid','on')
%     for k = 1:Nz
%         text(k,ad(k,1)+0.15*sign(ad(k,1)), ...
%             [num2str(nm(k,1)) ' ' num2str(nm(k,2))], ...
%             'HorizontalAlignment','center')
%     end
%     title(['Modified Zernike Spectrum (N = ' num2str(N) ')'])
%     movegui(h,'center'), set(h,'Visible','on')
%
%     % Display the (30) modified Zernike basis functions
%     z(is_in) = 0;
%     h = figure('Visible','off','Position',[0 0 1000 800]);
%     K = reshape(1:30,6,5)'; K = K(:);
%     for k = 1:30
%         zk = z;
%         zk(is_in) = zernfun(nm(k,1),nm(k,2),r,theta + ad(k,2));
%         subplot(5,6,K(k))
%         pcolor(x,x,zk), shading interp
%         set(gca,'XTick',[],'YTick',[])
%         axis square
%         text(-1,1.2,['[ ' num2str(nm(k,:)) ' ]'])
%         text(0.45,1.2,sprintf('%0.2f',ad(k,1)))
%     end
%     uicontrol('Style','text','Position',[0 770 1000 20], ...
%               'String','Modified Zernike Modal Contributions', ...
%               'HorizontalAlignment','center','FontSize',12)
%     movegui(h,'center'), set(h,'Visible','on')
%
%
%   Description
%   -----------
%   The standard Zernike basis is defined as
%
%      Znm(r,theta) = Rnm(r) * cos( m*theta)  m > 0
%      Zn0(r,theta) = Rn0(r)                           m = [-n:2:n]
%      Znm(r,theta) = Rnm(r) * sin(-m*theta)  m < 0
%
%   Any function F(r,theta) on a circular domain can be written as a sum
%   over this basis (i.e., over n and m),
%
%      F(r,theta) = sum{ Cnm * Znm(r,theta) }  .
%
%   For some systems, it is more convenient or informative to characterize
%   functions or data using the alternative (modified) Zernike basis 
%   (Campbell [1])
%
%      Ynm(r,theta) = Rnm(r) * cos(m*(theta - Dnm))  m = [0:2:n] for n even
%      Yn0(r,theta) = Rn0(r)                         m = [1:2:n] for n odd
%
%   for which
%
%      F(r,theta) = sum{ Anm * Ynm(r,theta;Dnm) }
%                 = sum{ Anm * Rnm(r) * cos(m*(theta - Dnm)) }  .
%
%   With the standard Zernike basis, each radial mode Rnm(r) or wavefront
%   aberration (such as astigmatism) is essentially decomposed into a pair
%   of mutually orthogonal components,
%
%      Cnm(+) * Rnm(r) * cos( m*theta)     m > 0
%      Cnm(-) * Rnm(r) * sin(-m*theta)     m < 0
% 
%   With the modififed basis, each mode is instead decribed as a single
%   term, with a single amplitude and a rotation angle in the (X,Y)-plane,
%
%      Anm * Rnm(r) * cos(m*(theta - Dnm))  .
%
%   The new modal coefficients are computed from those of the standard fit
%   as
%    
%      Anm = sqrt( [Cnm(+)]^2 + [Cnm(-)]^2 )
%      Dnm = atan( [Cnm(-)] / [Cnm(+)] )
% 
%   Note that, since each (+m)/(-m)-coefficient pair is combined to compute
%   a single Anm, the modified Zernike spectrum consists of [n m]-modes
%   with positive m only:
%
%      m = [0:2:n] for n even
%      m = [1:2:n] for n odd
%
%   Reference:
%   [1] A New Method for Describing the Aberrations of the Eye Using
%        Zernike Polynomials. C. E. Campbell. Optometry and Vision 
%        Science. 80(1), 79�83 (2003).
%
%   See also ZERNFUN

%   Paul Fricker 2/28/2012

% Check and prepare the inputs:
% -----------------------------
if ( ~any(size(r    )) == 1 ) || ...
   ( ~any(size(theta)) == 1 ) || ...
   ( ~any(size(data )) == 1 )
    error('zernmodfit:NMvectors1', ...
          'The inputs R, THETA, and DATA must all be vectors.')
end

if ( numel(r)     ~= numel(theta)) || ...
   ( numel(r)     ~= numel(data) ) || ...
   ( numel(theta) ~= numel(data) )
    error('zernmodfit:NMvectors2', ...
          'The inputs R, THETA, and DATA must all have the same number of elements.')
end

if numel(N) > 1
    error('zernmodfit:NMvectors4', ...
          'N must be a single number (scalar - positive integer or zero).')
end

if N < 0 || N ~= round(N)
    error('zernmodfit:NMvectors3', ...
          'N must be a positive integer or zero.')
end

r     = r(:);
theta = theta(:);
data  = data(:);

if any( r>1 | r<0 )
    error('zernmodfit:Rlessthan1','All R must be between 0 and 1.')
end

% Build vectors of mode numbers:
% % ------------------------------
% original code
% n = cellfun(@(x)(x*ones(1,x+1)),num2cell(0:N),'UniformOutput',false);
% n = [n{:}]';
% m = cellfun(@(x)([fliplr(-x:2:-1) fliplr(x:-2:0)]),num2cell(0:N),'UniformOutput',false);
% m = [m{:}]';

% modified code
n = cellfun(@(x)(x*ones(1,x+1)),num2cell(0:N),'UniformOutput',false);
n = [n{:}]';
m = cellfun(@(x)([(-x:2:-1) fliplr(x:-2:0)]),num2cell(0:N),'UniformOutput',false);
m = [m{:}]';


% Compute the required Zernike functions:
% ---------------------------------------
r = r(:);
theta = theta(:);
z = zernfun(n,m,r,theta);

% Compute the standard Zernike coefficients:
% ------------------------------------------
c = z\data(:);

% Assign the outputs as standard form:
% -------------------
ad = [c zeros(size(c))];
varargout = {[n m]};

% % Separate the positive, negative, and zero cases:
% % ------------------------------------------------
% cp = c(m>0);
% cn = c(m<0);
% cz = c(m==0);
%         
% % Compute the magnitudes and rotation angles from the coefficients:
% % -----------------------------------------------------------------
% cnz = sqrt(cp.^2 + cn.^2);
% dnz = atan2(-cn,cp);
% dnz(dnz<0) = dnz(dnz<0) + 2*pi;
% dnz = dnz./m(m>0);
% 
% % Combine the m==0 and m~=0 coefficients and angles into vectors:
% % ---------------------------------------------------------------
% a = zeros(size([cz;cnz]));
% d = a;
% m_pos_and_zero = m(m>=0);
% is_zero = m_pos_and_zero==0;
% is_pos  = m_pos_and_zero>0;
% a(is_zero) = cz;
% a(is_pos)  = cnz;
% d(is_zero) = 0;
% d(is_pos)  = dnz;
% 
% % Assign the outputs:
% % -------------------
% ad = [a d];
% varargout = {[n(m>=0) m(m>=0)]};