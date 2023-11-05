// Intercept only model

data {
	int<lower=1> N;
	vector[N] y; 		// response
	}

parameters {
  //real mu_a;   
  //real<lower=0> sigma_a;
  real<lower=0> sigma_y; 
  real a; // intercept
	}

transformed parameters {
   real yhat[N];
       	for(i in 1:N){
            yhat[i] = a;
		// b[i] * x[i];
			     	}

	}

model {
        a ~ normal(25, 2);
        sigma_y ~ normal(2, 2);
	//a ~ normal(mu_a, sigma_a); 	
        //mu_a ~ normal(0, 50);
        //sigma_a ~ normal(0, 10);
	
	y ~ normal(yhat, sigma_y);

}
