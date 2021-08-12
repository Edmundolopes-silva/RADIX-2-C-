#include <cassert>
#include <complex>
#include <vector>
// To demonstrate runtime
#include <chrono>
#include <iostream>

#define PI 3.1416

static std::vector<std::complex<double>> FFTRadix2(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& W);
static bool IsPowerOf2(size_t x);
static size_t ReverseBits(const size_t x, const size_t n);

std::vector<std::complex<double>> FFT(const std::vector<double>& x) //FUNÇÃO FFT
{
    size_t N = x.size();

    // Radix2 FFT exige que N seja uma potência de 2
    
    assert(IsPowerOf2(N)); 

    // por conta da simetria, dividimos em pares e impares e conseguimos computar a parte real com um FFT complexa de metade dos pontos, basta carregar depois da divisão em pares e impares para um vetor complexo
    std::vector<std::complex<double>> x_p(N / 2); // vetor x_p de (N/2) pontos 
    //std::cout<< "Vetor complexo:" <<std::endl;
    for (size_t n = 0; n < N / 2; ++n)
    {
        x_p[n] = std::complex<double>(x[2 * n], x[2 * n + 1]);
        //std:: cout<<  x_p[n] << std::endl;
    }

    // Fatores de torcao
    std::vector<std::complex<double>> W(N / 2); //FATORES IMPARES
    std::vector<std::complex<double>> W_p(N / 4);//FATORES PARES
    //std::cout<< "Fatores de torcao:" <<std::endl;
    for (size_t k = 0; k < N / 2; ++k)
    {
        W[k] = std::polar(1.0, -2 * M_PI * k / N);
       // std::cout<< W[k] <<std::endl;
        // O DFT complexo de N / 2 pontos usa apenas os fatores de torção pares
        if (k % 2 == 0)
        {
            W_p[k / 2] = W[k];
        }
    }
    //Fatores pares
    //std::cout<< "Fatores de torcao pares:" <<std::endl;
    for(size_t k = 0; k < N / 4; ++k)
    {
        //std::cout<< W_p[k] <<std::endl;
    }
    
   
    std::vector<std::complex<double>> X_p = FFTRadix2(x_p, W_p); //VETOR COMPLEXO DE N/2 PONTOS E FATORES PARES

    // Extraia o FFT do ponto N do sinal real dos resultados 
    std::vector<std::complex<double>> X(N);
    X[0] = X_p[0].real() + X_p[0].imag();
    for (size_t k = 1; k < N / 2; ++k)
    {
        // Extraia o FFT dos componentes pares
        auto A = std::complex<double>(
            (X_p[k].real() + X_p[N / 2 - k].real()) / 2,
            (X_p[k].imag() - X_p[N / 2 - k].imag()) / 2);

        // Extraia o FFT dos componentes impares
        auto B = std::complex<double>(
            (X_p[N / 2 - k].imag() + X_p[k].imag()) / 2,
            (X_p[N / 2 - k].real() - X_p[k].real()) / 2);

        // Soma dos resultados 
        X[k] = A + W[k] * B;
        X[k + N / 2] = A - W[k] * B;
    }

    return X;
}

static std::vector<std::complex<double>> FFTRadix2(const std::vector<std::complex<double>>& x, const std::vector<std::complex<double>>& W)
{
    size_t N = x.size();

    // Radix2 FFT exige que N seja uma potência de 2
    assert(IsPowerOf2(N));

    // quantos estágios o FFT deve calcular
    size_t stages = static_cast<size_t>(log2(N));

    // vetor de saída com índices invertidos de bits
    std::vector<std::complex<double>> X(N);
    for (size_t n = 0; n < N; ++n)
    {
        X[n] = x[ReverseBits(n, stages)];
    }

    // Calcule a FFT um estágio de cada vez e soma os resultados
    for (size_t stage = 1; stage <= stages; ++stage)
    {
        size_t N_stage = static_cast<size_t>(std::pow(2, stage));
        size_t W_offset = static_cast<size_t>(std::pow(2, stages - stage));
        for (size_t k = 0; k < N; k += N_stage)
        {
            for (size_t n = 0; n < N_stage / 2; ++n)
            {
                auto tmp = X[k + n];
                X[k + n] = tmp + W[n * W_offset] * X[k + n + N_stage / 2];
                X[k + n + N_stage / 2] = tmp - W[n * W_offset] * X[k + n + N_stage / 2];
            }
        }
    }

    return X;
}

// RETORNA TRUE SE FOR DE POTENCIA 2
static bool IsPowerOf2(size_t x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
    //return x && (!(x & (x - 1)));
}

// x é composto de n bits, reversão dos bits
static size_t ReverseBits(const size_t x, const size_t n)
{
    size_t xReversed = 0;
    for (size_t i = 0; i < n; i++)
    {
        xReversed = (xReversed << 1U) | ((x >> i) & 1U);
    }

    return xReversed;
}

int main()
{
    size_t N = 16;  //quantidade de pontos
    std::vector<double> x(N); //vetor x[n]
    std::vector<int> b(N);

    int f_s = 8000;  //FREQUENCIA de amostragem
    double t_s = 1.0 / f_s; // TAXA DE amostragem

    for (size_t n = 0; n < N; ++n)
    {
        x[n] = sin(2*PI*1000*n*t_s) + 0.5*sin(2*PI*2000*n*t_s + (3*PI/4)); //sinal
        
    }

    auto start = std::chrono::high_resolution_clock::now(); //INICIO DO TEMPO
    auto X = FFT(x); //RESULTADO DA FFT
    auto stop = std::chrono::high_resolution_clock::now(); //FINAL DO TEMPO
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start); //TEMPO EM microseconds
    for(int a=0; a<N; a++)
    {
       std::cout<< "PONTO " << a << " : " << X[a] << std::endl;
       
    }
    std::cout << "microseconds: " << duration.count() << std::endl;
}


