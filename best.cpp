#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <random>
#include <fstream>
#include <complex>

#include <chrono>

#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/special_functions/digamma.hpp>

//global variables
int M, iterations, n_chains, sample_freq, na_muji, na_N, na_phij, ploidy;
double mu, m_N, v_N, m_phij, v_phij, lnL, tunning_muji, tunning_N, tunning_phij;
std::vector<int> counts;


// g++ mcmc_neutral.cpp -o mcmc -O3
// ./mcmc sfs.txt 0.00015 100000 2 10 

// general functions
void             ReadControlFile           (const std::string str_control_file,std::string& model, std::string& str_sfile);
std::vector<int> s_Simulator               (double& muij, double& muji, double& phij, int& N, int& S);
void             s_SimulatorOddM           (double& muij, double& muji, double& phij, int& N, std::vector<double>& msfs, const std::vector<double>& nsfs);
void             s_SimulatorEvenM          (double& muij, double& muji, double& phij, int& N, std::vector<double>& msfs, const std::vector<double>& nsfs);

// functions for the Moran model with mutation biases and selection
double           s_Likelihood              (double& muij, double& muji, double& phij, int& N);
double           s_logLikelihoodOddM       (int& N, const std::vector<double>& nsfs, double& sum_nsfs);
double           s_logLikelihoodEvenM      (int& N, const std::vector<double>& nsfs, double& sum_nsfs);
void             s_MetropolisHastings_muji (double& muij, double& muji, double& phij, int& N);
void             s_MetropolisHastings_phij (double& muij, double& muji, double& phij, int& N);
void             s_MetropolisHastings_N    (double& muij, double& muji, double& phij, int& N);
void             s_MCMC                    (std::string& str_ifile);

// functions for the neutral Moran model with mutation biases
double           n_Likelihood              (double& muij, double& muji, int& N);
void             n_MetropolisHastings_muji (double& muij, double& muji, int& N);
void             n_MetropolisHastings_N    (double& muij, double& muji, int& N);
void             n_MCMC                    (std::string& str_ifile);

// statistical functions
double runif     (double min, double max);
double rnorm     (double mean, double sd);
double rexp      (double lambda);
double lbinomial (int n, int k);
double dhyper    (int v, int r, int n, int N);
double rgamma(double alpha, double beta);



int main(int argc, char * argv[]) {

  // gets the control file name  
  std::string str_control_file = argv[1]; 

  // reads control file and parameters related to the data, prior and mcmc
  std::string str_sfile;
  std::string model; 
  std::vector<int> counts;
  ReadControlFile(str_control_file,model,str_sfile);
  

  
  // simUlates a sample SFS
  int N    = 10000;
  double muij = 0.0002;
  double muji = 0.0001;
  double phij = 1.001;
  mu   = muji*muij*(1.0+pow(phij,N-1))/(muji+muij*pow(phij,N-1)) ;
  M       = 200;
  int S   = 1000000;

  counts = s_Simulator(muij,muji,phij,N,S);

  std::ofstream sfile(str_sfile);
  for (int i=0; i<(M+1); ++i){
    sfile << counts[i] << " ";
  }
  sfile.close();
  double lll = s_Likelihood(muij,muji,phij,N);
  std::cout << "lnL:" << lll << "\n";



  ReadControlFile(str_control_file,model,str_sfile);

  mu   = muji*muij*(1.0+pow(phij,N-1))/(muji+muij*pow(phij,N-1)) ;
  
  

  //mu = m_mu;
  //std::cout << "mu1:" << mu << "\n";

  // runs the neutral moran model
  if (model=="N") {
    std::cout << "  Parameters were succefully read.\n\n";
    std::cout << model << " | " << iterations << " | " << n_chains << " | " << sample_freq << " | "  << mu*10000 << " | " << m_N << " | " << v_N << " | " << m_phij << " | " << v_phij << "\n\n";
    n_MCMC(str_sfile);
    std::cout << "\n  The log and chekpoint files have been written.\n\n";

  // runs the moran model with selection
  } else if (model == "S") {
    std::cout << "  Parameters were succefully read.\n\n";
    std::cout << model << " | " << iterations << " | " << n_chains << " | " << sample_freq << " | "  << mu*10000 << " | " << m_N << " | " << v_N << " | " << m_phij << " | " << v_phij << "\n\n";
    s_MCMC(str_sfile);
    std::cout << "\n  The log and chekpoint files have been written.\n\n";

  } else {
    std::cerr << "  Model is incorrectly specified. Check control file.\n\n";
  }
  
  
}





/* 
  ReadControlFile
  
  Read the control file which specifies the counts file,
  the parameters of the parameters priors and the parameters to set MCMC

*/

void ReadControlFile(const std::string str_control_file, std::string& model, std::string& str_sfile) {

  //reads the control file 
  std::string word;

  std::ifstream cfile(str_control_file);
  if (!cfile.is_open()) {
    std::cerr << "  Unable to open " << str_control_file << ".\n\n";
  }

  // updates parameters
  cfile >> word >> model >> word >> str_sfile >> word >> mu >> word >> m_N >> word >> v_N >> word >> m_phij >> word >> v_phij >> word >> iterations >> word >> n_chains >> word >> sample_freq; 
  cfile.close();

  // getting the denominators of normal priors
  v_N    = 2*pow(v_N,2);
  v_phij = 2*pow(v_phij,2);

  // reads the sample site frequency spectrum
  std::ifstream sfile(str_sfile);
  if (!sfile.is_open()) {
    std::cerr << "  Unable to open " << str_sfile << ".\n\n";
  }

  // gets the empirical site frequency spectrum and builds counts vector
  int value;  
  while ( sfile >> value ) {
    counts.push_back(value);
  }
  sfile.close();

  // gets the sample size
  M = counts.size()-1;

}




/*
  Likelihood

  Calculates the likelihood of observing a vector of site 
  counts sampled from M individuals given that these follow a 
  Moran dynamic with N individuals and two alleles i and j 
  governed by boundary mutations (muji and muij) and 
  selection (phij)

  Arguments:
  muij    mutation rate from allele i to j
  muji    mutation rate from allele j to i
  phij    fitness of allele j (phii = 1)
  N       number of haploid invidials in the population
  M       number of sampled individuals
  counts  site counts observed in a sample of M individuals 
          following the order: {0j},{1j},...,{M-1j},{Mj}

  Output:
  The ln likelihood

*/

double s_Likelihood_Approximation(double& muij, double& muji, double& phij, int& N){

  int K=10000;
  double f,choose;
  
  double B[K+1][M+1];
  for (int m=0; m<(M+1); m++){
    choose = exp(lbinomial(M,m));
    for (int k=0; k<(K+1); ++k){
      f = k*1.0/K;
      B[k][m]= choose*pow(f,m)*pow(f,1.0*M-m);
      //B[k,m]= choose*pow(f,m)*pow(f,1.0*M-m);
    }
  }

  std::vector<double> nsfs(K+1);
  double sum_nsfs = 0.0;
  nsfs[0] = muji;
  nsfs[K] = muij*pow(phij,N-1);
  sum_nsfs += nsfs[0] + nsfs[K];
  for (int k=1; k<K; ++k){
    f = k*1.0/K;
    nsfs[k] = muij*muji*pow(phij,f*N-1)*(f*phij+1-f)/(f*(1-f)) ;
    sum_nsfs += nsfs[k];  
  }

  std::vector<double> esfs(M+1,0.0);
  double lnL=0.0;
  for (int m=0; m<(M+1); m++){
    for (int k=0; k<(K+1); ++k){
      esfs[m] += B[k][m]*nsfs[k];
    }
    lnL += counts[m]*(log(esfs[m])-log(sum_nsfs));
  }



}

double s_Likelihood(double& muij, double& muji, double& phij, int& N) {
  
  // computing the site frequency spectrum for a population of 
  // N individuals
  // {0j},{1j},...,{nj},...,{N-1j},{Nj}
  std::vector<double> nsfs(N+1);

  // observed fixed sites {0j} and {Mj}
  nsfs[0] = muji;
  nsfs[N] = muij*pow(phij,N-1);
  double sum_nsfs = nsfs[0] + nsfs[N];

  // polymorphic states {nj}
  double mus = N*muij*muji;
  double phij_0 = 1.0/phij;

  for (int n = 1; n<N; ++n){
    phij_0 *= phij;
    nsfs[n] = mus*phij_0*(n*phij+N-n)/(n*(N-n));
    sum_nsfs += nsfs[n];
  }

  //std::cout << "\nsum_nsfs:" << sum_nsfs << "\n";
  //for (int i=0; i<(N+1); ++i){
  //    std::cout << "nsfs_"<< i << ": " << nsfs[i] << "\n";
  //}

  // computing the sampled site frequency spectrum for a sample of 
  // M individuals from a population of N individuals
  // {0j},{1j},...,{mj},...,{M-1j},{Mj}
  if ( M%2==0 ) { 
    return s_logLikelihoodEvenM(N,nsfs,sum_nsfs); 
  } else { 
    return s_logLikelihoodOddM(N,nsfs,sum_nsfs); 
  }


}

double s_logLikelihoodOddM(int& N, const std::vector<double>& nsfs, double& sum_nsfs){

  // preliminaries: M is the sample size and N is the population size
  std::vector<double> msfs(M+1,0.0);
  int scount = 0;
  double hyper, base;
  double lnLikelihood = 0.0;
  double lbNM = lbinomial(N,M);

  // going through {0j} ... {mj} ... {Mj}
  for (int m=0; m<(M/2+1); ++m){
    
    base = lbinomial(N-m,M-m) - lbNM; 
    hyper = exp(base);
    msfs[m]   += nsfs[m]*hyper;
    msfs[M-m] += nsfs[N-m]*hyper; 
 
    // going through {0j} ... {nj} ... {Nj}
    for (int n=(m+1); n<(N-M+m+1); ++n){

      base += log(n) + log(N-n-M+m+1) - log(n-m) - log(N-n+1);
      hyper = exp(base);
      msfs[m]   += nsfs[n]*hyper;
      msfs[M-m] += nsfs[N-n]*hyper;
    }

    lnLikelihood += counts[m]*log(msfs[m])+counts[M-m]*log(msfs[M-m]);
    scount += counts[m] + counts[M-m];
  }

  lnLikelihood -= scount*log(sum_nsfs);

  return lnLikelihood;

  //for (int i=0; i<(M+1); ++i){
  //    std::cout << "msfs_"<< i << ": " << msfs[i] << "\n";
  //}
}

double s_logLikelihoodEvenM(int& N, const std::vector<double>& nsfs, double& sum_nsfs) {

  // preliminaries: M is the sample size and N is the population size
  std::vector<double> msfs(M+1,0.0);
  int scount = 0;
  double hyper, base;
  double lnLikelihood = 0.0;
  double lbNM = lbinomial(N,M);

  // going through {0j} ... {mj} ... {Mj}
  for (int m=0; m<(M/2); ++m){
    
    base = lbinomial(N-m,M-m) - lbNM; 
    hyper = exp(base);
    msfs[m]   += nsfs[m]*hyper;
    msfs[M-m] += nsfs[N-m]*hyper; 
 
    for (int n=(m+1); n<(N-M+m+1); ++n){

      base += log(n) + log(N-n-M+m+1) - log(n-m) - log(N-n+1);
      hyper = exp(base);
      msfs[m]   += nsfs[n]*hyper;
      msfs[M-m] += nsfs[N-n]*hyper;
    }

    lnLikelihood += counts[m]*log(msfs[m])+counts[M-m]*log(msfs[M-m]);
    scount += counts[m] + counts[M-m];
  }  

  base = lbinomial(N-M/2,M-M/2) - lbNM; 
  hyper = exp(base);
  msfs[M/2]   += nsfs[M/2]*hyper;
 
  for (int n=(M/2+1); n<(N-M+M/2+1); ++n){
    base += log(n) + log(N-n-M+M/2+1) - log(n-M/2) - log(N-n+1);
    hyper = exp(base);
    msfs[M/2]   += nsfs[n]*hyper;
  }

  scount += counts[M/2];
  lnLikelihood += counts[M/2]*log(msfs[M/2]) - scount*log(sum_nsfs);


  return lnLikelihood;

  //for (int i=0; i<(M+1); ++i){
  //    std::cout << "msfs_"<< i << ": " << msfs[i] << "\n";
  //}
}



/*
  MetropolisHastings_muji
  
  Performs a Metropolis-Hastings step on the mutation rate muji
  by employing the Metropolis-Hastings algorithm with a 
  Normal proposal

  Arguments:
  mu           average mutation rate = 0.5*muji+0.5*muij 
  lnL          initial likelihood
  tunning_muji controls the standard deviation of the normal 
               proposal
  na_muji      number of times a newly proposed muji was accepted 
               during the Metropolis-Hastings step

  Output:
  Updates muji, muij, lnL and na_muji internally
*/

void s_MetropolisHastings_muji(double& muij, double& muji, double& phij, int& N){
  
  // Normal proposal
  double muji1 = rnorm(muji,tunning_muji);
  double muij1 = muji1*mu/( muji1*(1.0+pow(phij,N-1))-pow(phij,N-1)*mu );  
  //std::cout << "muji: " << muji << " / miji1: " << muji1 << " : " << tunning_muji << "," << na_muji << "\n";

  if ( (muji1<0.0) || (muji1>1.0) ) {
    return; 
  }
  
  // Likelihood
  double lnL1 = s_Likelihood(muij1,muji1,phij,N);
  
  // Accept-reject step
  // proposal ratio = 1 
  double acceptance_prob = exp(lnL1-lnL);
  if (runif(0,1) < acceptance_prob) {
    muji     = muji1;
    muij     = muij1;
    lnL      = lnL1; 
    na_muji += 1;    
  } else {
    return; 
  }

}


/*
  MetropolisHastings_phij
  
  Performs a Metropolis-Hastings step on the fitness phij 
  by employing the Metropolis-Hastings algorithm with a 
  Normal proposal

  Arguments:
  lnL          initial likelihood
  tunning_phij controls the standard deviation of the normal 
               proposal
  na_phij      number of times a newly proposed phij was accepted 
               during the Metropolis-Hastings step

  Output:
  Updates phij, lnL and na_phij internally

*/

void s_MetropolisHastings_phij(double& muij, double& muji, double& phij, int& N){
  
  // Normal proposal
  double phij1 = rnorm(phij,tunning_phij);
  //std::cout << "phij: " << phij << " / phij1: " << phij1 << " : " << tunning_phij << "," << na_phij << "\n";

  if ( (phij1)<0.0 || ( (phij1-1.0) >1.0/N) ) {
    return; 
  }
  
  // Likelihood
  double lnL1 = s_Likelihood(muij,muji,phij1,N);
  
  // Accept-reject step
  // proposal ratio = 1 
  double acceptance_prob = exp( lnL1 -lnL -pow(phij1-m_phij,2)/v_phij +pow(phij-m_phij,2)/v_phij );
  if (runif(0,1) < acceptance_prob) {
    phij     = phij1;
    lnL      = lnL1; 
    na_phij += 1;    
  } else {
    return; 
  }

}





/*
  MetropolisHastings_Nij
  
  Performs a Metropolis-Hastings step on the effective population
  size N by employing the Metropolis-Hastings algorithm with a 
  Normal proposal

  Arguments:
  lnL          initial likelihood
  tunning_N    controls the standard deviation of the normal 
               proposal
  na_N         number of times a newly proposed N was accepted 
               during the Metropolis-Hastings step

  Output:
  Updates N, lnL and na_N internally

*/

void s_MetropolisHastings_N(double& muij, double& muji, double& phij, int& N){
  
  // Normal proposal
  int N1 = rnorm(N*1.0,tunning_N);
  //std::cout << "N: " << N << " / N1: " << N1 << " : " << tunning_N << "," << na_N << "\n";

  if ( (N1<(M+1)) || (N1==N) ) {
    return; 
  }
  
  // Likelihood
  double lnL1 = s_Likelihood(muij,muji,phij,N1);

  // Accept-reject step
  // proposal ratio = 1 
  double acceptance_prob = exp( lnL1 -lnL -pow(N1-m_N,2)/v_N +pow(N-m_N,2)/v_N );
  if (runif(0,1) < acceptance_prob) {
    N     = N1;
    lnL   = lnL1; 
    na_N += 1;    
  } else {
    return; 
  }

}


/*
  MCMC

  Performs a MCMC 


*/

void s_MCMC(std::string& str_ifile){

  std::string str_ofile,str_ckfile;

  // employing the mcmc independent n_chains
  for (int c=0; c<n_chains; ++c){

    // creating log and checkpoint output files
    str_ofile = str_ifile + std::to_string(c+1) + ".log";
    str_ckfile = str_ifile + std::to_string(c+1) + ".checkpoint";
    std::ofstream ofile(str_ofile);
    std::ofstream ckfile(str_ckfile);

    // preliminaries
    ofile << "gen\tlnL\tmuji\tmuji\tphij\tN\n"; 
    ckfile << "gen\tmuji\tmuij\tphij\tN\tna_muji\tn_phij\tna_N\ttunning_muji\ttunning_phij\ttunning_N\n";

    // sampling the initial parameters
    double muji = rnorm(mu,0.01*mu);
    int N       = round(rexp(mu) + counts.size()-1);
    //N         = 1000;
    double phij = 1.0+1.0/(2.0*N); // rexp(N);
    //phij = 1.00001;

    double muij = muji*mu/( muji*(1.0+pow(phij,N-1))-pow(phij,N-1)*mu );  
    //std::cout << "rmuji: " << muji << "  rmuij: " << muij << "\n"; 

    // setting the acceptance probabilitites and initial tunning
    // parameters
    na_muji = 0; 
    na_phij = 0; 
    na_N = 0;
    tunning_muji = mu;
    tunning_phij = 1.0/N;
    tunning_N    = N;

    // initial likelihood
    lnL = s_Likelihood(muij,muji,phij,N);

    // output
    std::cout << 0 << "\t" << lnL << "\t" << muji << "\t" << muij << "\t" << phij << "\t" << N << "\n";

    for (int i=1; i<(iterations+1); ++i){

      // Metropolis-Hastings algorithm
      s_MetropolisHastings_muji(muij,muji,phij,N);
      s_MetropolisHastings_phij(muij,muji,phij,N);
      s_MetropolisHastings_N(muij,muji,phij,N);
      //std::cout << i+1 << "\t" << lnL << "\t" << muji << "\t" << muij << "\t" << phij << "\t" << N << "\t" << muji+muij<< "\n";

      if (i%sample_freq == 0) {

        // tunning the proposal variances
        tunning_muji /= (2.0*0.234/(0.234+1.0*na_muji/i));
        tunning_phij /= (2.0*0.234/(0.234+1.0*na_phij/i));
        tunning_N    /= (2.0*0.234/(0.234+1.0*na_N/i));
      
        //std::cout << i+1 << "\t" << na_muji << "\t" << na_phij << "\t" << na_N << "\t" << tunning_muji << "\t" << tunning_phij << "\t" << tunning_N << "\n";

        std::ostringstream os_muij,os_muji,os_lnL,os_tunning_muji,os_tunning_phij; 
        os_muij << muij;
        os_muji << muji;
        os_lnL  << lnL;
        os_tunning_muji << tunning_muji;
        os_tunning_phij << tunning_phij;

        // saving the sampled parameters
        ofile << i << "\t" << os_lnL.str() << "\t" << os_muji.str() << "\t" << os_muij.str() << "\t" << std::to_string(phij) << "\t" << std::to_string(N) << "\n"; 
        ckfile << i << "\t" << os_muij.str() << "\t" << os_muji.str() << "\t" << std::to_string(phij) << std::to_string(N) << "\t" << std::to_string(na_muji) << "\t" <<  std::to_string(na_phij) << "\t" << std::to_string(na_N) << "\t" << os_tunning_muji.str() << "\t" << os_tunning_phij.str() << "\t" << std::to_string(tunning_N) << "\n"; 
      
        //output
        std::cout << i+1 << "\t" << lnL << "\t" << muji << "\t" << muij << "\t" << phij << "\t" << N << "\n";

      }
    }
  
    // closes
    ofile.close();
    ckfile.close();
  }

}


double n_Likelihood(double& muij, double& muji, int& N){

  // computing the sampled site frequency spectrum for a sample of 
  // M individuals from a population of N individuals
  // {0j},{1j},...,{mj},...,{M-1j},{Mj}
  std::vector<double> msfs(M+1,0.0);

  // observed fixed sites {0j} and {Mj}
  msfs[0] = muji+muij*muji*N*(boost::math::digamma(N) - boost::math::digamma(M));
  msfs[M] = muij+muij*muji*N*(boost::math::digamma(N) - boost::math::digamma(M));
  double sum_msfs = muji+muij*muji*N*(boost::math::digamma(N) - boost::math::digamma(M)) + muij+muij*muji*N*(boost::math::digamma(N) - boost::math::digamma(M));

  // observed polymorphic sites {mj}
  for (int m = 1; m<M; ++m){
    msfs[m] = muij*muji*N*M/(m*(M-m));
    sum_msfs += msfs[m]; 
  }

  // calculating the log likelihood 
  double lnLikelihood = 0.0;
  // expected number of sites per site pattern {mj}
  for (int m = 0; m<(M+1); ++m) {
    lnLikelihood +=  counts[m]*log(msfs[m]/sum_msfs);
  }

  return lnLikelihood;
}

void n_MetropolisHastings_muji(double& muij, double& muji, int& N){

  // Normal proposal
  double muji1 = rnorm(muji,tunning_muji);
  double muij1 = mu*muji1/(2.0*muji1-mu);

  if ( (muji1<0.0) || (muji1>1.0) ) {
    return; 
  }
  
  // Likelihood
  double lnL1 = n_Likelihood(muij1,muji1,N);
  //std::cout << "lnl_mu: " << lnL1 << "\n";

  // Accept-reject step
  // proposal ratio = 1 
  double acceptance_prob = exp(lnL1-lnL);
  if (runif(0,1) < acceptance_prob) {
    muji     = muji1;
    muij     = muij1;
    lnL      = lnL1; 
    na_muji += 1;    
  } else {
    return; 
  }

}

void n_MetropolisHastings_N(double& muij, double& muji, int& N){

  // Normal proposal
  int N1 = rnorm(N*1.0,tunning_N);
  
  if ( (N1<(M+1)) || (N1==N) ) {
    //std::cout << "Reject: " << N << " / "<< N1 << "\n";
    return; 
  }
  
  // Likelihood
  double lnL1 = n_Likelihood(muij,muji,N1);
  //std::cout << "lnl_N: " << lnL1 << "\n";

  // Accept-reject step
  // proposal ratio = 1 
  double acceptance_prob = exp( lnL1 -lnL -pow(N1-m_N,2)/v_N +pow(N-m_N,2)/v_N );
  //std::cout << "\nap: "<< acceptance_prob << "\n" ;
  if (runif(0,1) < acceptance_prob) {
    N     = N1;
    lnL   = lnL1; 
    na_N += 1;    
  } else {
    return; 
  }
}


void n_MCMC(std::string& str_ifile){

  std::string str_ofile,str_ckfile;

  // employing the mcmc independent n_chains
  for (int c=0; c<n_chains; ++c){

    // creating log and checkpoint output files
    str_ofile = str_ifile + std::to_string(c+1) + ".log";
    str_ckfile = str_ifile + std::to_string(c+1) + ".checkpoint";
    std::ofstream ofile(str_ofile);
    std::ofstream ckfile(str_ckfile);

    // preliminaries
    ofile << "gen\tlnL\tmuji\tmuji\tN\n"; 
    ckfile << "gen\tmuji\tmuij\tN\tna_muji\tna_N\ttunning_muji\ttunning_N\n";

    // sampling the initial parameters
    double muji = rnorm(mu,mu*0.01);
    int N       = round(rexp(0.001) + counts.size()-1);
    //N       = 1000;

    double muij = mu*muji/(2.0*muji-mu);

    // setting the acceptance probabilitites and initial tunning
    // parameters
    na_muji = 0; 
    na_N = 0;
    tunning_muji = mu;
    tunning_N    = N;

    // initial likelihood
    lnL = n_Likelihood(muij,muji,N);

    // output
    std::cout << 0 << "\t" << lnL << "\t" << muji << "\t" << muij << "\t" << N << "\n";

    for (int i=1; i<(iterations+1); ++i){

      // Metropolis-Hastings algorithm
      n_MetropolisHastings_muji(muij,muji,N);
      n_MetropolisHastings_N(muij,muji,N);
      //std::cout << i+1 << "\t" << lnL << "\t" << muji << "\t" << muij << "\t" << phij << "\t" << N << "\t" << muji+muij<< "\n";

      if (i%sample_freq == 0) {

        // tunning the proposal variances
        tunning_muji /= (2.0*0.234/(0.234+1.0*na_muji/i));
        tunning_N    /= (2.0*0.234/(0.234+1.0*na_N/i));
      
        //std::cout << i+1 << "\t" << na_muji << "\t" << na_phij << "\t" << na_N << "\t" << tunning_muji << "\t" << tunning_phij << "\t" << tunning_N << "\n";

        std::ostringstream os_muij,os_muji,os_lnL,os_tunning_muji; 
        os_muij << muij;
        os_muji << muji;
        os_lnL  << lnL;
        os_tunning_muji << tunning_muji;

        // saving the sampled parameters
        ofile  << i << "\t" << os_lnL.str()  << "\t" << os_muji.str() << "\t" << os_muij.str()     << "\t" << std::to_string(N)       << "\n"; 
        ckfile << i << "\t" << os_muij.str() << "\t" << os_muji.str() << "\t" << std::to_string(N) << "\t" << std::to_string(na_muji) << "\t" << std::to_string(na_N) << "\t" << os_tunning_muji.str() << "\t" << std::to_string(tunning_N) << "\n"; 
      
        //output
        std::cout << i+1 << "\t" << lnL << "\t" << muji << "\t" << muij << "\t" << N << "\n";

      }
    }
  
    // closes
    ofile.close();
    ckfile.close();
  }

}



/*
  Simulator

  Simulates the site frequency spectrum of a sample of M 
  individuals and S sites from a Moran population with two 
  alleles i and j and N haploid individuals governed by 
  boundary mutations (muij and muji) and selection (phij). 

  Arguments:
  muij  mutation rate from allele i to j
  muji  mutation rate from allele j to i
  phij  fitness of allele j (phii = 1)
  N     number of haploid invidials in the population
  M     number of sampled individuals
  S     number of sampled sites
  
  Output:
  A vector of sampled site counts following the order:
  {0j},{1j},...,{mj},...,{M-1j},{Mj}

*/

std::vector<int> s_Simulator(double& muij, double& muji, double& phij, int& N, int& S) {

  // computing the site frequency spectrum for a population of 
  // N individuals
  // {0j},{1j},...,{nj},...,{N-1j},{Nj}
  std::vector<double> nsfs(N+1);

  // observed fixed sites {0j} and {Mj}
  nsfs[0] = muji;
  nsfs[N] = muij*pow(phij,N-1);
  double sum_nsfs = nsfs[0] + nsfs[N];

  // polymorphic states {nj}
  double mus = N*muij*muji;
  double phij_0 = 1.0/phij;

  for (int n = 1; n<N; ++n){
    phij_0 *= phij;
    nsfs[n] = mus*phij_0*(n*phij+N-n)/(n*(N-n));
    sum_nsfs += nsfs[n];
  }


  std::vector<double> msfs(M+1,0.0);

  if ( M%2==0 ) { s_SimulatorEvenM(muij,muji,phij,N,msfs,nsfs); }
  else { s_SimulatorOddM(muij,muji,phij,N,msfs,nsfs); }

  // expected number of sites per site pattern {mj}
  std::vector<int> ssfs(M+1,0.0);
  double factor = S/sum_nsfs;
  for (int m = 0; m<(M+1); ++m) {
    ssfs[m] = round(msfs[m]*factor);
    std::cout << ssfs[m] << "\n";
  }

  return ssfs;

}


void s_SimulatorOddM(double& muij, double& muji, double& phij, int& N,std::vector<double>& msfs, const std::vector<double>& nsfs){

  // preliminaries: M is the sample size and N is the population size
  int scount = 0;
  double hyper, base;
  double lbNM = lbinomial(N,M);

  // going through {0j} ... {mj} ... {Mj}
  for (int m=0; m<(M/2+1); ++m){
    
    base = lbinomial(N-m,M-m) - lbNM; 
    hyper = exp(base);
    msfs[m]   += nsfs[m]*hyper;
    msfs[M-m] += nsfs[N-m]*hyper; 
 
    // going through {0j} ... {nj} ... {Nj}
    for (int n=(m+1); n<(N-M+m+1); ++n){

      base += log(n) + log(N-n-M+m+1) - log(n-m) - log(N-n+1);
      hyper = exp(base);
      msfs[m]   += nsfs[n]*hyper;
      msfs[M-m] += nsfs[N-n]*hyper;
    }
  }

}


void s_SimulatorEvenM(double& muij, double& muji, double& phij, int& N,std::vector<double>& msfs, const std::vector<double>& nsfs) {

  // preliminaries: M is the sample size and N is the population size
  int scount = 0;
  double hyper, base;
  double lbNM = lbinomial(N,M);

  // going through {0j} ... {mj} ... {Mj}
  for (int m=0; m<(M/2); ++m){
    
    base = lbinomial(N-m,M-m) - lbNM; 
    hyper = exp(base);
    msfs[m]   += nsfs[m]*hyper;
    msfs[M-m] += nsfs[N-m]*hyper; 
 
    for (int n=(m+1); n<(N-M+m+1); ++n){

      base += log(n) + log(N-n-M+m+1) - log(n-m) - log(N-n+1);
      hyper = exp(base);
      msfs[m]   += nsfs[n]*hyper;
      msfs[M-m] += nsfs[N-n]*hyper;
    }

  }  

  base = lbinomial(N-M/2,M-M/2) - lbNM; 
  hyper = exp(base);
  msfs[M/2]   += nsfs[M/2]*hyper;
 
  for (int n=(M/2+1); n<(N-M+M/2+1); ++n){
    base += log(n) + log(N-n-M+M/2+1) - log(n-M/2) - log(N-n+1);
    hyper = exp(base);
    msfs[M/2]   += nsfs[n]*hyper;
  }


}



// Other usefull statistical functions


double runif(double min, double max){
  std::random_device rd; 
  std::mt19937 generator(rd()); 
  std::uniform_real_distribution<double> distribution(min,max);
  double random = distribution(generator);
  return random;
}

double rnorm(double mean, double sd){
  std::random_device rd; 
  std::mt19937 generator(rd()); 
  std::normal_distribution<double> distribution(mean,sd);
  double random = distribution(generator);
  return random;
}

double rexp(double lambda){
  std::random_device rd; 
  std::mt19937 generator(rd()); 
  std::exponential_distribution<double> distribution(lambda);
  double random = distribution(generator);
  return random;
}

double dhyper(int v, int r, int n, int N){
  boost::math::hypergeometric hyper(r, n, N);
  double probability = pdf(hyper, v);
  return probability;
}

double rgamma(double alpha, double beta){
  std::random_device rd; 
  std::mt19937 generator(rd()); 
  std::gamma_distribution<double> distribution(alpha,beta);
  double random = distribution(generator);
  return random;
}

double lbinomial(int n, int k){

  if (k == 0) { return 1; }
  else if (2*k > n) { lbinomial(n,n-k); }
  else {
    double e = log(n-k);
    for (int i=2; i<(k+1); ++i) {
      e += log(n-k+i) -log(i);
    }
    return e;
  }
}

