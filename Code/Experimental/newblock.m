function newblock(n,p)
	[K0,K1]=generateKsm(p);
	K=generateMatrix(n,K0,K1);
	[vK eK]=eig(K,'vector');
	[M0,M1]=generateMsm(p);
	M=generateMatrix(n,M0,M1);
	[vM eM]=eig(M,'vector');
	eL=eig(K,M);

	% symbol
	fK=@(t)K0+(K1+K1')*cos(t)+1i*((K1-K1')*sin(t));
	fM=@(t)M0+(M1+M1')*cos(t)+1i*((M1-M1')*sin(t));
	fL=@(t)fM(t)\fK(t);
	t=linspace(0,pi,n+1)';
	aeKs=[];aeMs=[];aeLs=[];
	for jj=1:n+1
		aeKs=[aeKs sort(eig(fK(t(jj))))];
		aeMs=[aeMs sort(eig(fM(t(jj))))];
		aeLs=[aeLs sort(eig(fL(t(jj))))];
	end

	% K
	aeK{1}=aeKs(1,2:n)';
	if p>1
		for qq=2:p
			bq=mod(qq,2);
			aeK{qq}=aeKs(qq,2-bq:n+1-bq)';
		end
	end

	% M
	aeM{1}=aeMs(1,2:n)';
	if p>1
		b=mod(p,2);
		if b==0 || p==3 
			for qq=2:p
				bq=mod(qq,2);
				aeM{qq}=aeMs(qq,2-bq:n+1-bq)';  % good
				% aeM{qq}=aeMs(qq,2:n+1)';  % not correct
			end
		else
			phat=p-2+2*mod((p+1)/2,2);
			for qq=2:(phat+1)/2
				bq=mod(qq,2);
				aeM{qq}=aeMs(qq,1+bq:n+bq)';
			end
			for qq=(phat+1)/2+1:p
				bq=mod(qq,2);
				aeM{qq}=aeMs(qq,2-bq:n+1-bq)';
			end
		end
	end

	%L
	b=mod(p,2);
	aeL{1}=aeLs(1,2:n+1-b)';
	if p>1
		for qq=2:2:p
			aeL{qq}=aeLs(qq,2-b:n+b)';
		end
		for qq=3:2:p
			aeL{qq}=aeLs(qq,1+b:n+1-b)';
		end
	end

	saeK=[];saeM=[];saeL=[];
	for qq=1:p
		saeK=[saeK;aeK{qq}];
		saeM=[saeM;aeM{qq}];
		saeL=real([saeL;aeL{qq}]);
	end
	[eK saeK]
	[eM saeM]
	[eL saeL]
	norm(eM-sort(saeM))
% 	size(vK)
% % 	figure, hold on 
% % 	for ii=30:39
% % 		plot(vK(:,ii))
% % end
% 	% for ii=1:n
% 	% 	plot(vK(:,ii+n),'--')
% 	% end
% 	N=n*p-1;
% 	tt=linspace(pi/(N+1),n*pi/(N+1),N)';
% 	eid=5;
% 	for kk=1:N
% 		av(kk)=sin(kk*tt(2*eid+2));
% 	end
% 	av=av'/norm(av);
% 	[vK(:,eid) av]
% figure, hold on 
% plot(vK(:,eid))
% plot(av)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=generateMatrix(n,T0,T1)
	% generates a tridiagonal block matrix
	T=kron(eye(n),T0)+...
	  kron(diag(ones(1,n-1),-1),T1)+...
	  kron(diag(ones(1,n-1), 1),T1');
	T=T(1:end-1,1:end-1);
end

function [K0,K1]=generateKsm(p)
	% K sub matrices for order p
	K0=zeros(p);
	K1=zeros(p);
	for ii=1:p
		for jj=1:p
			f=@(t)mydL(ii,t,p).*mydL(jj,t,p);
			K0(ii,jj)=integral(f,0,1);
		end
	end
	f=@(t)mydL(0,t,p).*mydL(0,t,p);
	K0(p,p)=K0(p,p)+integral(f,0,1);
	for ii=1:p
		f=@(t)mydL(0,t,p).*mydL(ii,t,p);
		K1(ii,p)=integral(f,0,1);
	end
end

function [M0,M1]=generateMsm(p)
	% M submatrices for order p
	M0=zeros(p);
	M1=zeros(p);
	for ii=1:p
		for jj=1:p
			f=@(t)myL(ii,t,p).*myL(jj,t,p);
			M0(ii,jj)=integral(f,0,1);
		end
	end
	f=@(t)myL(0,t,p).*myL(0,t,p);
	M0(p,p)=M0(p,p)+integral(f,0,1);
	for ii=1:p
		f=@(t)myL(0,t,p).*myL(ii,t,p);
		M1(ii,p)=integral(f,0,1);
	end
end


function L=myL(h,t,p)
	L=1;
	for k=0:p
		if k~=h
			L=L.*(t-k/p)./(h/p-k/p);
		end
	end
end
function dL=mydL(h,t,p)
	dL=0;
	for i=0:p
		if i~=h
			tmp=1;
			for k=0:p
				if k~=h && k~=i
					tmp=tmp.*(t-k/p)./(h/p-k/p);
				end
			end
			dL=dL+1/(h/p-i/p)*tmp;
		end
	end
end