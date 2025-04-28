function [X,P,meanInX,meanInY,varInX,varInY] = elliptic_cylinder_update(X,P,Y,sig_r,lambda,wm,wc,nxu,n_upd,meanInX,meanInY,varInX,varInY,tao,artificial_noise,filter,source)
    % dimensions
    nx = length(X);
    nMeas = size(Y,2);
    
    % check if number of measurements smaller than n_upd
    alpha = 0.1; beta = 2; kappa = 0;
    if nMeas < n_upd
        n_upd = nMeas;
        if strcmp(filter,'ERHM')
            nxu = nx + n_upd;
        elseif strcmp(filter,'GAM')
            nxu = nx;
        end
        lambda = alpha^2*(nxu + kappa) - nxu;
        wm = zeros(1,2*nxu+1); wc = wm;
        wm(1) = lambda/(nxu + lambda);
        wc(1) = lambda/(nxu + lambda) + (1 - alpha^2 + beta);
        wm(2:2*nxu + 1) = 1/(2*(nxu + lambda));
        wc(2:2*nxu + 1) = 1/(2*(nxu + lambda));
    end
    
    % iterate n times over n_upd measurements
    for j = 1:ceil(nMeas/n_upd)
        % get new measurement set
        if j <= floor(nMeas/n_upd)
            meas = Y(:,n_upd*j-n_upd+1:n_upd*j);
        else
            meas = Y(:,n_upd*(j-1)+1:end);
            % change number of measurements to be updated and UKF weights
            n_upd = mod(nMeas,n_upd);
            if strcmp(filter,'ERHM')
                nxu = nx + n_upd;
            elseif strcmp(filter,'GAM')
                nxu = nx;
            end
            lambda = alpha^2*(nxu + kappa) - nxu;
            wm = zeros(1,2*nxu+1); wc = wm;
            wm(1) = lambda/(nxu + lambda);
            wc(1) = lambda/(nxu + lambda) + (1 - alpha^2 + beta);
            wm(2:2*nxu + 1) = 1/(2*(nxu + lambda));
            wc(2:2*nxu + 1) = 1/(2*(nxu + lambda));
        end

        % calculate sigma points
        if strcmp(filter,'ERHM')
            Xu = [X; zeros(n_upd,1)];
            L = blkdiag(chol(P)',eye(n_upd));
        elseif strcmp(filter,'GAM')
            Xu = X;
            L = chol(P)';
        else
            error('Wrong filter setting!')
        end
        A = sqrt(nxu + lambda) * L;
        sig_points = [zeros(size(Xu)) -A A];
        sig_points = sig_points + repmat(Xu, 1, size(sig_points, 2));
        sig_points(nx+1:end,:) = 0.5*(1 + erf(sig_points(nx+1:end,:)/sqrt(2)));
        numSamples = size(sig_points,2);
    
        % predict sigma points
        z_predict = zeros(3*n_upd,numSamples);
        for i = 1:numSamples
            % extract sigma point parameters
            pos = sig_points(1:3,i);
            or = sig_points(4,i);
            a = c1(sig_points(5,i),0,'lower');
            b = c1(sig_points(6,i),0,'lower');
            h = c1(sig_points(7,i),0,'lower');
            % transform measurements to local coordinate system
            R_rot = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
            z_loc = R_rot'*(meas - pos);
            z_loc(3,:) = abs(z_loc(3,:));
            % extrusion factors
            if strcmp(filter,'ERHM')
                us = sig_points(8:end,i);
            elseif strcmp(filter,'GAM')
                us = min(max(2*z_loc(3,:)/h,0),1)';
            else
                us = 0;
            end
            if strcmp(source,'radial')
                % get angle parameter of measurement source
                theta = atan2(z_loc(2,:)*a,z_loc(1,:)*b);
                % measurement prediction
                zp = [a*cos(theta); b*sin(theta); us'*h/2];
            elseif strcmp(source,'projected')
                % polygonal chain
                ps = as_polygon(a,b,25);
                % project measurement to ellipse
                zp = [project(ps, z_loc(1:2,:)); us'*h/2];
            else
                error('Wrong source setting!')
            end
            z_predict(:,i) = zp(:) - z_loc(:);
        end
    
        % get predicted measurements
        z_pred = sum(wm.*z_predict,2);
    
        % calculate measurement noise covariance matrix
        if artificial_noise
            % for asymmetric noise
            % first step measurements in local coordinates
            pos = X(1:3); or = X(4); 
            a = c1(X(5),0,'lower'); b = c1(X(6),0,'lower'); 
            R_rot = [cos(or) -sin(or) 0; sin(or) cos(or) 0; 0 0 1];
            meas_loc = R_rot'*(meas - pos);
            % second step measurements inside or outside
            inside_outside = meas_loc(1,:).^2/(a - 3*sig_r)^2 + meas_loc(2,:).^2/(b - 3*sig_r)^2 - 1;
            % indices of measurements inside
            idx = find(inside_outside < 0);  
            % update estimate of inside measurement mean
            meanInX = 1/(1 + length(idx)/tao)*meanInX + 1/(tao + length(idx))*sum(abs(z_pred(3*idx-2)));
            meanInY = 1/(1 + length(idx)/tao)*meanInY + 1/(tao + length(idx))*sum(abs(z_pred(3*idx-1)));
            z_pred(3*idx-2) = z_pred(3*idx-2) + meanInX; z_pred(3*idx-1) = z_pred(3*idx-1) + meanInY;
            % update estimate of inside measurement variance
            varInX = 1/(1 + length(idx)/tao)*varInX + 1/(tao + length(idx))*sum((z_pred(3*idx-2) - meanInX).^2);
            varInY = 1/(1 + length(idx)/tao)*varInY + 1/(tao + length(idx))*sum((z_pred(3*idx-1) - meanInY).^2);
            % build covariance matrix
            R = sig_r^2*ones(1,3*n_upd);
            R(3*idx-2) = varInX; R(3*idx-1) = varInY; R = diag(R);
        else
            R = sig_r^2*eye(3*n_upd);
        end
    
        % get update matrices
        S = zeros(3*n_upd);
        Psi = zeros(nx,3*n_upd);
        for i = 1:size(sig_points,2)
            S = S + wc(i)*(z_predict(:,i) - z_pred)*(z_predict(:,i) - z_pred)';
            Psi = Psi + wc(i)*(sig_points(1:nx,i) - X)*(z_predict(:,i) - z_pred)';
        end
        S = S + R;
        
        % do UKF measurement update
        K = Psi/S;
        X = X + K*(zeros(3*n_upd,1) - z_pred);
        P = P - K*Psi';
        P = 0.5*(P+P');
    end
end