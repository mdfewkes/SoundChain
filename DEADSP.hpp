#pragma once

#include <memory>
#include <cmath>
#include <cstring>

namespace DEA {

	class DSP_Bypass {
	public:
		struct Parameters {};

		DSP_Bypass() {}
		DSP_Bypass(Parameters &parameters) : _params(parameters ) {}
		~DSP_Bypass() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;
		}

		void Reset(int sampleRate) {}

		float GetSample(float xn) {
			return xn;
		}

	private:
		Parameters _params;
	};

	class DSP_Biquad {
	public:
		enum FilterType {LowPass, HighPass, 
		                 BandPass, BandStop, 
		                 LowShelf, HighShelf,
		                 ButterworthLowPass, ButterworthHighPass,
		                 LinkwitzRileyLowPass, LinkwitzRileyHighPass,
		                 AllPass1stOrder, AllPass2ndOrder,
		                 LowPass1stOrder, HighPass1stOrder,
		                 Peak};
		struct Parameters {
			float frequency = 1000.0f;
			float q = 0.707f;
			float filterGain = 0.0f;
			FilterType filterType = FilterType::LowPass;
		};
		struct Coefficients {
			float a0 = 1.0f;
			float a1 = 0.0f;
			float a2 = 0.0f;
			float b1 = 0.0f;
			float b2 = 0.0f;
			float c0 = 1.0f;
			float d0 = 0.0f;
		};

		DSP_Biquad() {}
		DSP_Biquad(Parameters &parameters) : _params(parameters ) {}
		~DSP_Biquad() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;
			CalculateCoefficients();
		}

		Coefficients GetCoefficients() {return _coeffs;}
		void SetCoefficients(Coefficients &coefficients) {
			_coeffs = coefficients;
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;

			_x_z1 = 0.0f;
			_x_z2 = 0.0f;
			_y_z1 = 0.0f;
			_y_z2 = 0.0f;

			CalculateCoefficients();
		}

		float GetSample(float xn) {
			float yn = _coeffs.a0 * xn
			         + _coeffs.a1 * _x_z1 
			         + _coeffs.a2 * _x_z2 
			         - _coeffs.b1 * _y_z1
			         - _coeffs.b2 * _y_z2;

			_x_z2 = _x_z1;
			_x_z1 = xn;

			_y_z2 = _y_z1;
			_y_z1 = yn;

			return yn * _coeffs.c0 + xn * _coeffs.d0;
		}

	private:
		Parameters _params;
		Coefficients _coeffs;
		int _sampleRate = -1;
		float _x_z1 = 0.0f;
		float _x_z2 = 0.0f;
		float _y_z1 = 0.0f;
		float _y_z2 = 0.0f;

		void CalculateCoefficients() {
			switch (_params.filterType) {
				case FilterType::LowPass:
					CalculateLowPassCoefficients();
					break;
				case FilterType::HighPass:
					CalculateHighPassCoefficients();
					break;
				case FilterType::BandPass:
					CalculateBandPassCoefficients();
					break;
				case FilterType::BandStop:
					CalculateBandStopCoefficients();
					break;
				case FilterType::LowShelf:
					CalculateLowShelfCoefficients();
					break;
				case FilterType::HighShelf:
					CalculateHighShelfCoefficients();
					break;
				case FilterType::ButterworthLowPass:
					CalculateButterworthLowPassCoefficients();
					break;
				case FilterType::ButterworthHighPass:
					CalculateButterworthHighPassCoefficients();
					break;
				case FilterType::LinkwitzRileyLowPass:
					CalculateLinkwitzRileyLowPassCoefficients();
					break;
				case FilterType::LinkwitzRileyHighPass:
					CalculateLinkwitzRileyHighPassCoefficients();
					break;
				case FilterType::AllPass1stOrder:
					CalculateAllPass1stOrderCoefficients();
					break;
				case FilterType::AllPass2ndOrder:
					CalculateAllPass2ndOrderCoefficients();
					break;
				case FilterType::LowPass1stOrder:
					CalculateLowPass1stOrderCoefficients();
					break;
				case FilterType::HighPass1stOrder:
					CalculateHighPass1stOrderCoefficients();
					break;
				case FilterType::Peak:
					CalculatePeakCoefficients();
					break;
			} 
		}

		void CalculateLowPassCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / _sampleRate;
			float cosw0 = cos(w0);
			float sinw0 = sin(w0);
			float d = 1.0f / _params.q;
			float beta = 0.5f * ((1.0f - d / 2.0f * sinw0) / (1.0f + d / 2.0f * sinw0));
			float gamma = (0.5f + beta) * cosw0;

			_coeffs.a0 = (0.5f + beta - gamma) / 2.0f;
			_coeffs.a1 = 0.5f + beta - gamma;
			_coeffs.a2 = (0.5f + beta - gamma) / 2.0f;
			_coeffs.b1 = -2.0f * gamma;
			_coeffs.b2 = 2.0f * beta;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateHighPassCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / _sampleRate;
			float cosw0 = cos(w0);
			float sinw0 = sin(w0);
			float d = 1.0f / _params.q;
			float beta = 0.5f * ((1.0f - d / 2.0f * sinw0) / (1.0f + d / 2.0f * sinw0));
			float gamma = (0.5f + beta) * cosw0;

			_coeffs.a0 = (0.5f + beta + gamma) / 2.0f;
			_coeffs.a1 = -(0.5f + beta + gamma);
			_coeffs.a2 = (0.5f + beta + gamma) / 2.0f;
			_coeffs.b1 = -2.0f * gamma;
			_coeffs.b2 = 2.0f * beta;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateBandPassCoefficients() {
			float k = tan((M_PI * _params.frequency) / _sampleRate);
			float k2 = k*k;
			float delta = k2 * _params.q + k + _params.q;

			_coeffs.a0 = k / _params.q;
			_coeffs.a1 = 0.0f;
			_coeffs.a2 = -k / delta;
			_coeffs.b1 = (2.0f * _params.q * (k2 - 1)) / delta;
			_coeffs.b2 = (k2 * _params.q - k + _params.q) / delta;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateBandStopCoefficients() {
			float k = tan((M_PI * _params.frequency) / _sampleRate);
			float k2 = k*k;
			float delta = k2 * _params.q + k + _params.q;

			_coeffs.a0 = (_params.q * (k2 + 1)) / delta;
			_coeffs.a1 = (2.0f * _params.q * (k2 - 1)) / delta;
			_coeffs.a2 = (_params.q * (k2 + 1)) / delta;
			_coeffs.b1 = (2.0f * _params.q * (k2 - 1)) / delta;
			_coeffs.b2 = (k2 * _params.q - k + _params.q) / delta;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateLowShelfCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / _sampleRate;
			float mu = std::pow(10, _params.filterGain / 20.0f);
			float beta = 4.0f / (1.0f + mu);
			float delta = beta * tan(w0 / 2.0f);
			float gamma = (1.0f - delta) / (1.0f + delta);

			_coeffs.a0 = (1 - gamma) / 2.0f;
			_coeffs.a1 = (1 - gamma) / 2.0f;
			_coeffs.a2 = 0.0f;
			_coeffs.b1 = -gamma;
			_coeffs.b2 = 0.0f;
			_coeffs.c0 = mu - 1.0f;
			_coeffs.d0 = 1.0f;
		}

		void CalculateHighShelfCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / _sampleRate;
			float mu = std::pow(10, _params.filterGain / 20.0f);
			float beta = (1.0f + mu) / 4.0f;
			float delta = beta * tan(w0 / 2.0f);
			float gamma = (1.0f - delta) / (1.0f + delta);

			_coeffs.a0 = (1 + gamma) / 2.0f;
			_coeffs.a1 = (1 + gamma) / -2.0f;
			_coeffs.a2 = 0.0f;
			_coeffs.b1 = -gamma;
			_coeffs.b2 = 0.0f;
			_coeffs.c0 = mu - 1.0f;
			_coeffs.d0 = 1.0f;
		}

		void CalculateButterworthLowPassCoefficients() {
			float c = 1.0f / tan((M_PI * _params.frequency) / _sampleRate);
			float a0 = 1.0f / (1.0f + std::sqrt(2.0f * c) + c*c);

			_coeffs.a0 = a0;
			_coeffs.a1 = 2.0f * a0;
			_coeffs.a2 = a0;
			_coeffs.b1 = 2.0f * a0 * (1.0f - c*c);
			_coeffs.b2 = a0 * (1.0f - std::sqrt(2.0f * c) + c*c);
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateButterworthHighPassCoefficients() {
			float c = tan((M_PI * _params.frequency) / _sampleRate);
			float a0 = 1.0f / (1.0f + std::sqrt(2.0f * c) + c*c);

			_coeffs.a0 = a0;
			_coeffs.a1 = -2.0f * a0;
			_coeffs.a2 = a0;
			_coeffs.b1 = 2.0f * a0 * (c*c - 1.0f);
			_coeffs.b2 = a0 * (1.0f - std::sqrt(2.0f * c) + c*c);
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateLinkwitzRileyLowPassCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / _sampleRate;
			float omega = M_PI * _params.frequency;
			float kappa = omega / tan(w0);
			float delta = kappa*kappa + omega*omega + 2.0f * kappa * omega;
			float omega2 = omega*omega;
			float kappa2 = kappa*kappa;

			_coeffs.a0 = omega2 / delta;
			_coeffs.a1 = 2.0f * (omega2 / delta);
			_coeffs.a2 = omega2 / delta;
			_coeffs.b1 = (-2.0f * kappa2 + 2.0f * omega2) / delta;
			_coeffs.b2 = (-2.0f * kappa * omega + kappa2 + omega2) / delta;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateLinkwitzRileyHighPassCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / _sampleRate;
			float omega = M_PI * _params.frequency;
			float kappa = omega / tan(w0);
			float delta = kappa*kappa + omega*omega + 2.0f * kappa * omega;
			float omega2 = omega*omega;
			float kappa2 = kappa*kappa;

			_coeffs.a0 = kappa2 / delta;
			_coeffs.a1 = (-2.0f * kappa2) / delta;
			_coeffs.a2 = kappa2 / delta;
			_coeffs.b1 = (-2.0f * kappa2 + 2.0f * omega2) / delta;
			_coeffs.b2 = (-2.0f * kappa * omega + kappa2 + omega2) / delta;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateAllPass1stOrderCoefficients() {
			float w0 = M_PI * _params.frequency / _sampleRate;
			float tanw0 = tan(w0);
			float alpha = (tanw0 - 1) / (tanw0 +1);

			_coeffs.a0 = alpha;
			_coeffs.a1 = 1.0f;
			_coeffs.a2 = 0.0f;
			_coeffs.b1 = alpha;
			_coeffs.b2 = 0.0f;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateAllPass2ndOrderCoefficients() {
			float bw = _params.frequency / _params.q;
			float tanbwpisr = tan(bw * M_PI / _sampleRate);
			float alpha = (tanbwpisr - 1) / (tanbwpisr +1);
			float beta = -cos(2.0f * M_PI * _params.frequency / _sampleRate);

			_coeffs.a0 = -alpha;
			_coeffs.a1 = beta * (1.0f - alpha);
			_coeffs.a2 = 1.0f;
			_coeffs.b1 = beta * (1.0f - alpha);
			_coeffs.b2 = -alpha;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateLowPass1stOrderCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / _sampleRate;
			float cosw0 = cos(w0);
			float sinw0 = sin(w0);
			float gamma = cosw0 / (1.0f + sinw0);

			_coeffs.a0 = (1.0f - gamma) / 2.0f;
			_coeffs.a1 = (1.0f - gamma) / 2.0f;
			_coeffs.a2 = 0.0f;
			_coeffs.b1 = -gamma;
			_coeffs.b2 = 0.0f;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculateHighPass1stOrderCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / _sampleRate;
			float cosw0 = cos(w0);
			float sinw0 = sin(w0);
			float gamma = cosw0 / (1.0f + sinw0);

			_coeffs.a0 = (1.0f + gamma) / 2.0f;
			_coeffs.a1 = -((1.0f + gamma) / 2.0f);
			_coeffs.a2 = 0.0f;
			_coeffs.b1 = -gamma;
			_coeffs.b2 = 0.0f;
			_coeffs.c0 = 1.0f;
			_coeffs.d0 = 0.0f;
		}

		void CalculatePeakCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / _sampleRate;
			float mu = std::pow(10, _params.filterGain / 20.0f);
			float zeta = 4.0f / (1.0f + mu);
			float zetatanw02q = zeta * tan(w0 / 2.0f * _params.q);
			float beta = 0.5f * (1.0f - zetatanw02q) / (1.0f + zetatanw02q);
			float gamma = (0.5f + beta) * cos(w0);

			_coeffs.a0 = 0.5f - beta;
			_coeffs.a1 = 0.0f;
			_coeffs.a2 = -(0.5f - beta);
			_coeffs.b1 = -2.0f * gamma;
			_coeffs.b2 =2.0f * beta;
			_coeffs.c0 = mu - 1.0f;
			_coeffs.d0 = 1.0f;
		}
	};

	class DSP_RegaliaMitraFilter {
	public:
		enum FilterOrder {First, Second};
		struct Parameters {
			float frequency = 500.0f;
			float k = 0.65f;
			FilterOrder filterOrder = FilterOrder::Second;
		};

		DSP_RegaliaMitraFilter() {}
		DSP_RegaliaMitraFilter(Parameters &parameters) : _params(parameters ) {}
		~DSP_RegaliaMitraFilter() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			auto afpParams = _apf.GetParameters();
			afpParams.frequency = _params.frequency;
			afpParams.q = 1.0f;
			if (_params.filterOrder == FilterOrder::First) {
				afpParams.filterType = DSP_Biquad::FilterType::AllPass1stOrder;
			} else {
				afpParams.filterType = DSP_Biquad::FilterType::AllPass2ndOrder;
			}
			_apf.SetParameters(afpParams);
		}

		void Reset(int sampleRate) {
			_apf.Reset(sampleRate);

			SetParameters(_params);
		}

		float GetSample(float xn) {
			float apf = _apf.GetSample(xn);
			float yn = (xn + apf) * 0.5f + (xn - apf) * _params.k/2.0f;

			return yn;
		}

	private:
		Parameters _params;
		DSP_Biquad _apf;
	};

	class DSP_MutiFilter {
	public:
		enum FilterOrder {First, Second};
		enum OutputMode {LowPass, AllPass, HighPass};
		struct Parameters {
			float frequency = 1000.0f;
			float q = 0.707f;
			FilterOrder filterOrder = FilterOrder::Second;
			OutputMode outputMode = OutputMode::LowPass;
		};
		struct Output {
			float lpf;
			float apf;
			float hpf;
		};

		DSP_MutiFilter() {}
		DSP_MutiFilter(Parameters &parameters) : _params(parameters ) {}
		~DSP_MutiFilter() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			auto lfpParams = _lpf.GetParameters();
			lfpParams.frequency = _params.frequency;
			lfpParams.q = _params.q;
			if (_params.filterOrder == FilterOrder::First) {
				lfpParams.filterType = DSP_Biquad::FilterType::LowPass1stOrder;
			} else {
				lfpParams.filterType = DSP_Biquad::FilterType::LowPass;
			}
			_lpf.SetParameters(lfpParams);
		}

		Output GetOutput() {return _output;}

		void Reset(int sampleRate) {
			_lpf.Reset(sampleRate);
			SetParameters(_params);
		}

		float GetSample(float xn) {
			float yn;

			_output.lpf = _lpf.GetSample(xn);
			_output.hpf = xn - _output.lpf;
			_output.apf = _output.lpf - _output.hpf;

			switch (_params.outputMode) {
				case OutputMode::LowPass:
					yn = _output.lpf;
					break;
				case OutputMode::AllPass:
					yn = _output.apf;
					break;
				case OutputMode::HighPass:
					yn = _output.hpf;
					break;
			}

			return yn;
		}

	private:
		Parameters _params;
		DSP_Biquad _lpf;
		Output _output;
	};

	class DSP_LFO {
	public:
		enum WaveformType {Triangle, Sin, Saw};
		struct Parameters {
			float frequency = 1.0f;
			WaveformType waveformType = WaveformType::Sin;
		};
		struct Output {
			float normalOutput;
			float quadPhaseOutput;
			float invertedOutput;
			float negativeQuadPhaseOutput;
		};

		DSP_LFO() {}
		DSP_LFO(Parameters &parameters) : _params(parameters ) {}
		~DSP_LFO() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;
			_phaseInc = _params.frequency / (float)_sampleRate;

			_modCounter = 0.0f;
			_modCounterQP = 0.25f;
		}

		Output GetSample() {
			Output output;

			switch (_params.waveformType) {
				case WaveformType::Triangle:
					output.normalOutput = abs(_modCounter * 2.0f - 1.0f) * 2.0f - 1.0f;
					output.quadPhaseOutput = abs(_modCounterQP * 2.0f - 1.0f) * 2.0f - 1.0f;
					output.invertedOutput = -output.normalOutput;
					output.negativeQuadPhaseOutput = -output.quadPhaseOutput;
					break;
				case WaveformType::Sin:
					output.normalOutput = parabolicSine(-(_modCounter * 2.0f * M_PI - M_PI));
					output.quadPhaseOutput = parabolicSine(-(_modCounterQP * 2.0f * M_PI - M_PI));
					output.invertedOutput = -output.normalOutput;
					output.negativeQuadPhaseOutput = -output.quadPhaseOutput;
					break;
				case WaveformType::Saw:
					output.normalOutput = _modCounter * 2.0f - 1.0f;
					output.quadPhaseOutput = _modCounterQP * 2.0f - 1.0f;
					output.invertedOutput = -output.normalOutput;
					output.negativeQuadPhaseOutput = -output.quadPhaseOutput;
					break;
			}

			_modCounter += _phaseInc;
			if (_modCounter >= 1.0f) _modCounter -= 1.0f;
			_modCounterQP += _modCounter + 0.25f;
			if (_modCounterQP >= 1.0f) _modCounterQP -= 1.0f;

			return output;
		}

	private:
		Parameters _params;
		int _sampleRate = -1;
		float _modCounter;
		float _modCounterQP;
		float _phaseInc;

		float parabolicSine(float angle) {
			float B = 4.0 / M_PI;
			float C = -4.0 / (M_PI*M_PI);
			float P = 0.225;
			float x = angle;

			float y = B * x + C * x * abs(x);
			y = P * (y * abs(y) - y) + y;
			return y;
		}
	};

	class DSP_AudioDetector {
	public:
		enum AudioDetectMode {Peak, MS, RMS};
		struct Parameters {
			float attackTime_mSec = 0.0f;
			float releaseTime_mSec = 0.0f;
			AudioDetectMode detectMode = AudioDetectMode::Peak;
			bool detectDB = false;
			bool clampToUnityMax = true;
		};

		DSP_AudioDetector() {}
		DSP_AudioDetector(Parameters &parameters) : _params(parameters ) {}
		~DSP_AudioDetector() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			if (_params.attackTime_mSec != parameters.attackTime_mSec) {
				_attackTime = exp(ANAOG_TC / (parameters.attackTime_mSec * _sampleRate * 0.001f));
			}
			if (_params.releaseTime_mSec != parameters.releaseTime_mSec) {
				_releaseTime = exp(ANAOG_TC / (parameters.releaseTime_mSec * _sampleRate * 0.001f));
			}
			
			_params = parameters;
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;
			_lastValue = 0.0f;

			_attackTime = exp(ANAOG_TC / (_params.attackTime_mSec * _sampleRate * 0.001f));
			_releaseTime = exp(ANAOG_TC / (_params.releaseTime_mSec * _sampleRate * 0.001f));

			SetParameters(_params);
		}

		float GetSample(float xn) {
			float input = xn;

			// Rectifier
			input = std::abs(input);
			if (_params.detectMode != AudioDetectMode::Peak) { // MS or RMS
				input *= input;
			}

			// RC Simulation
			float currentValue = 0.0f;
			if (input > _lastValue) {
				currentValue = _attackTime * ( _lastValue - input) + input;
			} else {
				currentValue = _releaseTime * ( _lastValue - input) + input;
			}
			if (_params.clampToUnityMax) {
				currentValue = std::min(currentValue, 1.0f);
			}
			currentValue = std::max(currentValue, 0.0f);
			_lastValue = currentValue;

			// RMS and log conversion
			if (_params.detectMode == AudioDetectMode::RMS) {
				currentValue = std::pow(currentValue, 0.5f);
			}
			if (_params.detectDB) {
				if (currentValue <= 0.0f) {
					currentValue = -96.0f;
				} else {
					currentValue = 20.0f * std::log10(currentValue);
				}
			}

			return currentValue;
		}

	private:
		const float ANAOG_TC = -0.999672;

		Parameters _params;
		int _sampleRate = -1;
		float _attackTime = 0.0f;
		float _releaseTime = 0.0f;
		float _lastValue = 0.0f;
	};

	class DSP_EnvelopeFollower {
	public:
		struct Parameters {
			float frequency = 1000.0f;
			float q = 0.707f;
			float attackTime_mSec = 10.0f;
			float releaseTime_mSec = 10.0f;
			float threshold = 0.0f;
			float sensitivity = 1.0f;
		};

		DSP_EnvelopeFollower() {
			SetParameters(_params);
		}
		DSP_EnvelopeFollower(Parameters &parameters) : _params(parameters ) {
			SetParameters(_params);
		}
		~DSP_EnvelopeFollower() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			auto filterParams = _filter.GetParameters();
			filterParams.frequency = _params.frequency;
			filterParams.q = _params.q;
			filterParams.filterType = DSP_Biquad::FilterType::LowPass;
			_filter.SetParameters(filterParams);

			auto detectorParams = _detector.GetParameters();
			detectorParams.attackTime_mSec = _params.attackTime_mSec;
			detectorParams.releaseTime_mSec = _params.releaseTime_mSec;
			detectorParams.detectMode = DSP_AudioDetector::AudioDetectMode::RMS;
			detectorParams.detectDB = true;
			detectorParams.clampToUnityMax = false;
			_detector.SetParameters(detectorParams);

			_thresholdValue = pow(10.0f, _params.threshold / 20.0f);
		}

		void Reset(int sampleRate) {
			_filter.Reset(sampleRate);
			_detector.Reset(sampleRate);
			_maxFilterFrequency = sampleRate / 2.0f;
		}

		float GetSample(float xn) {
			float detect_dB = _detector.GetSample(xn);
			float detectValue = pow(10.0f, detect_dB / 20.0f);
			float delta = detectValue - _thresholdValue;

			if (delta > 0.0f) {
				auto filterParams = _filter.GetParameters();
				filterParams.frequency = _params.frequency;

				float modulatorValue = delta * _params.sensitivity;
				filterParams.frequency = modulatorValue * (_maxFilterFrequency - filterParams.frequency) + filterParams.frequency;

				_filter.SetParameters(filterParams);
			}

			return _filter.GetSample(xn);
		}

	private:
		Parameters _params;
		DSP_Biquad _filter;
		DSP_AudioDetector _detector;
		float _maxFilterFrequency = -1;
		float _thresholdValue = 0;
	};

	template <typename T>
	class CircularBuffer {
	public:
		CircularBuffer() {
			_buffer.reset(new T[0]);
		}
		~CircularBuffer() {}

		unsigned int GetBufferLength() {return _bufferLength;}

		void CreateBuffer(unsigned int bufferLength) {
			_writeIndex = 0;
			_bufferLength = (unsigned int)(pow(2, ceil(log(bufferLength)/log(2))));
			_wrapMask = _bufferLength - 1;
			_buffer.reset(new T[_bufferLength]);
			FlushBuffer();
		}

		void FlushBuffer() {
			for (int i = 0; i < _bufferLength; i++) {
				_buffer[i] = 0;
			}
		}

		void WriteBuffer(T input) {
			_buffer[_writeIndex++] = input;
			_writeIndex &= _wrapMask;
		}

		T ReadBuffer(int delayInSamples) {
			int readIndex = _writeIndex - delayInSamples;
			readIndex &= _wrapMask;
			return _buffer[readIndex];
		}

		T ReadBuffer(float delayInFractionalSamples, bool interpolate = true) {
			T start = ReadBuffer((int)delayInFractionalSamples);
			if (!interpolate) return start;

			T end = ReadBuffer((int)delayInFractionalSamples + 1);
			float t = delayInFractionalSamples - (int)delayInFractionalSamples;

			return Lerp(start, end, t);
		}

	private:
		std::unique_ptr<T[]> _buffer = nullptr;
		unsigned int _writeIndex;
		unsigned int _bufferLength;
		unsigned int _wrapMask;

		T Lerp (T start, T end, float t) {
			if (t >= 1.0f) return end;

			return t * end + (1.0f - t) * start;
		}
	};

	template <typename T>
	class LinearBuffer {
	public:
		LinearBuffer() {
			_buffer.reset(new T[0]);
		}
		~LinearBuffer() {}

		unsigned int GetBufferLength() {return _bufferLength;}

		void CreateBuffer(unsigned int bufferLength) {
			_bufferLength = (unsigned int)(pow(2, ceil(log(bufferLength)/log(2))));
			_buffer.reset(new T[_bufferLength]);
			FlushBuffer();
		}

		void FlushBuffer() {
			for (int i = 0; i < _bufferLength; i++) {
				_buffer[i] = 0;
			}
		}

		void WriteBuffer(size_t index, T input) {
			if (index >= _bufferLength) {
				index = index % _bufferLength;
			}

			_buffer[index] = input;
		}

		T ReadBuffer(size_t index) {
			if (index >= _bufferLength) {
				index = index % _bufferLength;
			}

			return _buffer[index];
		}

		T ReadBuffer(float index, bool interpolate = true) {
			if (index >= _bufferLength || index < 0.0f) return 0;

			T start = ReadBuffer((int)index);
			if (!interpolate) return start;

			T end = ReadBuffer((int)index + 1);
			float t = index - (int)index;

			return Lerp(start, end, t);
		}

	private:
		std::unique_ptr<T[]> _buffer = nullptr;
		unsigned int _bufferLength;

		T Lerp (T start, T end, float t) {
			if (t >= 1.0f) return end;

			return t * end + (1.0f - t) * start;
		}
	};

	class DSP_Delay {
	public:
		struct Parameters {
			float delayTime = 0.2f;
			float dry = 1.0f;
			float wet = 1.0f;
			float feedback = 0.5f;
		};

		DSP_Delay() {}
		DSP_Delay(Parameters &parameters) : _params(parameters ) {}
		~DSP_Delay() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			_delayInSamples = _params.delayTime * _sampleRate;
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;

			_delayInSamples = _params.delayTime * _sampleRate;

			_buffer.CreateBuffer(_sampleRate * 2);
		}

		float GetSample(float xn) {
			float yn = _buffer.ReadBuffer(_delayInSamples);
			float dn = xn + _params.feedback * yn;
			_buffer.WriteBuffer(dn);

			return _params.dry * xn + _params.wet * yn;
		}

	private:
		Parameters _params;
		CircularBuffer<float> _buffer;
		int _sampleRate;
		float _delayInSamples;
	};

	class DSP_AnalogDelay {
	public:
		struct Parameters {
			float delayTime = 0.2f;
			float dry = 1.0f;
			float wet = 1.0f;
			float feedback = 0.5f;
			float lpfFrequency = 1000.0f;
		};

		DSP_AnalogDelay() {}
		DSP_AnalogDelay(Parameters &parameters) : _params(parameters ) {}
		~DSP_AnalogDelay() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			_delayInSamples = _params.delayTime * _sampleRate;

			auto lpfParams = _lpf.GetParameters();
			lpfParams.frequency = _params.lpfFrequency;
			_lpf.SetParameters(lpfParams);
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;

			_delayInSamples = _params.delayTime * _sampleRate;

			_buffer.CreateBuffer(_sampleRate * 2);

			_lpf.Reset(sampleRate);
			auto lpfParams = _lpf.GetParameters();
			lpfParams.filterType = DSP_Biquad::FilterType::LowPass1stOrder;
			lpfParams.frequency = _params.lpfFrequency;
			_lpf.SetParameters(lpfParams);
		}

		float GetSample(float xn) {
			float yn = _lpf.GetSample(_buffer.ReadBuffer(_delayInSamples));
			float dn = xn + _params.feedback * yn;
			_buffer.WriteBuffer(dn);

			return _params.dry * xn + _params.wet * yn;
		}

	private:
		Parameters _params;
		CircularBuffer<float> _buffer;
		DSP_Biquad _lpf;

		int _sampleRate;
		float _delayInSamples;
	};

	class DSP_MultitapDelay {
	public:
		struct Parameters {
			float delayTimeTap1 = 0.15f;
			float delayTimeTap2 = 0.2f;
			float delayTimeTap3 = 0.3f;
			float delayTimeTap4 = 0.7f;
			float dry = 1.0f;
			float wet = 1.0f;
			float feedback = 0.5f;
			float lpfFrequency = 1000.0f;
		};

		DSP_MultitapDelay() {}
		DSP_MultitapDelay(Parameters &parameters) : _params(parameters) {}
		~DSP_MultitapDelay() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			_delayInSamples4 = _params.delayTimeTap1 * _sampleRate;
			_delayInSamples3 = _params.delayTimeTap2 * _sampleRate;
			_delayInSamples2 = _params.delayTimeTap3 * _sampleRate;
			_delayInSamples1 = _params.delayTimeTap4 * _sampleRate;

			auto lpfParams = _lpf.GetParameters();
			lpfParams.frequency = _params.lpfFrequency;
			_lpf.SetParameters(lpfParams);
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;

			_delayInSamples4 = _params.delayTimeTap1 * _sampleRate;
			_delayInSamples3 = _params.delayTimeTap2 * _sampleRate;
			_delayInSamples2 = _params.delayTimeTap3 * _sampleRate;
			_delayInSamples1 = _params.delayTimeTap4 * _sampleRate;

			_buffer.CreateBuffer(_sampleRate * 2);

			_lpf.Reset(sampleRate);
			auto lpfParams = _lpf.GetParameters();
			lpfParams.filterType = DSP_Biquad::FilterType::LowPass1stOrder;
			lpfParams.frequency = _params.lpfFrequency;
			_lpf.SetParameters(lpfParams);
		}

		float GetSample(float xn) {
			float mx = _buffer.ReadBuffer(_delayInSamples1)
			         + _buffer.ReadBuffer(_delayInSamples2)
			         + _buffer.ReadBuffer(_delayInSamples3)
			         + _buffer.ReadBuffer(_delayInSamples4);
			float yn = _lpf.GetSample(mx);
			float dn = xn + _params.feedback * yn;
			_buffer.WriteBuffer(dn);

			return _params.dry * xn + _params.wet * yn;
		}

	private:
		Parameters _params;
		CircularBuffer<float> _buffer;
		DSP_Biquad _lpf;

		int _sampleRate;
		float _delayInSamples1;
		float _delayInSamples2;
		float _delayInSamples3;
		float _delayInSamples4;
	};

	class DSP_ImpulseConvolver {
	public:
		static const unsigned int INIT_BUFFER_LENGTH = 512;

		struct Parameters {};

		DSP_ImpulseConvolver() {}
		DSP_ImpulseConvolver(Parameters &parameters) : _params(parameters ) {}
		~DSP_ImpulseConvolver() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;
		}

		void setImpulseResponce(float* irArray, unsigned int lengthPowerOfTwo) {
			if (lengthPowerOfTwo != _irLength) {
				_signalBuffer.CreateBuffer(lengthPowerOfTwo);
				_irBuffer.CreateBuffer(lengthPowerOfTwo);
				_irLength = lengthPowerOfTwo;
			}

			for (int i = 0; i < _irLength; i++) {
				_irBuffer.WriteBuffer((size_t)i, irArray[i]);
			}
		}

		unsigned int GetIRLength() {return _irLength;};
		unsigned int GetIRBufferLength() {return _irBuffer.GetBufferLength();};

		void Reset(int sampleRate) {
			if (_initialized) {
				_signalBuffer.FlushBuffer();
			} else {
				_signalBuffer.CreateBuffer(INIT_BUFFER_LENGTH);
				_irBuffer.CreateBuffer(INIT_BUFFER_LENGTH);
				_irLength = _irBuffer.GetBufferLength();
			}
		}

		float GetSample(float xn) {
			float yn = 0.0f;
			_signalBuffer.WriteBuffer(xn);

			// Convolution!
			float signal = 0.0f;
			float ir = 0.0f;
			for (int i = 0; i < _irLength; i++) {
				signal = _signalBuffer.ReadBuffer(i);
				ir = _irBuffer.ReadBuffer((size_t)i);
				yn += signal*ir;
			}

			return yn;
		}

	private:
		Parameters _params;
		CircularBuffer<float> _signalBuffer;
		LinearBuffer<float> _irBuffer;
		unsigned int _irLength = 0;
		bool _initialized = false;
	};

	class DSP_AnalogFIRFilter {
	public:
		enum FilterType {LPF1, HPF1, LPF2, HPF2, BPF2, BSF2};
		struct Parameters {
			FilterType filterType = FilterType::LPF2;
			float frequency = 1000.0f;
			float q = 0.707f;
		};

		DSP_AnalogFIRFilter() {}
		DSP_AnalogFIRFilter(Parameters &parameters) : _params(parameters ) {}
		~DSP_AnalogFIRFilter() {
			delete[] _magArray;
			delete[] _irArray;
		}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			calculateAnalogMagArray();
			freqSample(IR_LENGTH, _magArray, _irArray, true);
			_convolver.setImpulseResponce(_irArray, IR_LENGTH);
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;

			_convolver.Reset(sampleRate);
		}

		float GetSample(float xn) {
			return _convolver.GetSample(xn);
		}

	private:Parameters _params;
		const unsigned int IR_LENGTH = DSP_ImpulseConvolver::INIT_BUFFER_LENGTH;
		int _sampleRate = -1;
		float* _magArray = new float[IR_LENGTH];
		float* _irArray = new float[IR_LENGTH];
		unsigned int _dftArrayLen = IR_LENGTH;
		bool _mirrorMag = false;
		DSP_ImpulseConvolver _convolver;

		inline void calculateAnalogMagArray() {
			// --- calculate first half of array
			float actualLength = _mirrorMag ? (float)_dftArrayLen : (float)_dftArrayLen * 2.0f;
			uint32_t dumpLength = _mirrorMag ? _dftArrayLen / 2 : _dftArrayLen;

			float bin1 = _sampleRate / actualLength;// (float)_dftArrayLen;
			float zeta = 1.0f / (2.0f * _params.q);
			float w_c = 2.0f * M_PI * _params.frequency;

			// --- zero out array; if filter not supported, this will return a bank of 0's!
			memset(&_magArray[0], 0, _dftArrayLen * sizeof(float));

			for (uint32_t i = 0; i < dumpLength; i++) {
				float eval_w = 2.0f * M_PI * i * bin1;
				float w_o = eval_w / w_c;

				if (_params.filterType == FilterType::LPF1) {
					float denXSq = 1.0f + (w_o*w_o);
					_magArray[i] = 1.0f / (std::pow(denXSq, 0.5f));

				} else if (_params.filterType == FilterType::HPF1) {
					float denXSq = 1.0f + (w_o*w_o);
					_magArray[i] = w_o / (std::pow(denXSq, 0.5f));

				} else if (_params.filterType == FilterType::LPF2) {
					float denXSq = (1.0f - (w_o*w_o))*(1.0f - (w_o*w_o)) + 4.0f * zeta*zeta * w_o*w_o;
					_magArray[i] = 1.0f / (std::pow(denXSq, 0.5f));

				} else if (_params.filterType == FilterType::HPF2) {
					float denXSq = (1.0f - (w_o*w_o))*(1.0f - (w_o*w_o)) + 4.0f * zeta*zeta * w_o*w_o;
					_magArray[i] = (w_o*w_o) / (std::pow(denXSq, 0.5f));

				} else if (_params.filterType == FilterType::BPF2) {
					float denXSq = (1.0f - (w_o*w_o))*(1.0f - (w_o*w_o)) + 4.0f * zeta*zeta * w_o*w_o;
					_magArray[i] = 2.0f*w_o*zeta / (std::pow(denXSq, 0.5f));

				} else if (_params.filterType == FilterType::BSF2) {
					float numXSq = (1.0f - (w_o*w_o))*(1.0f - (w_o*w_o));
					float denXSq = (1.0f - (w_o*w_o))*(1.0f - (w_o*w_o)) + 4.0f * zeta*zeta * w_o*w_o;
					_magArray[i] = (std::pow(numXSq, 0.5f)) / (std::pow(denXSq, 0.5f));
				}
			}

			if (_mirrorMag) {
				int index = _dftArrayLen / 2 - 1;
				for (unsigned int i = _dftArrayLen / 2; i <_dftArrayLen; i++) {
					_magArray[i] = _magArray[index--];
				}
			}
		}

		inline void freqSample(int numCoeffs, float freqResponce[], float output[], int symmetry) {
			const float kTwoPi = 2.0f * M_PI;
			int n, k;
			float x, value, M;

			M = (numCoeffs - 1.0) / 2.0;
			if (symmetry) {
				if (numCoeffs % 2) {
					for (n = 0; n<numCoeffs; n++) {
						value = freqResponce[0];
						x = kTwoPi * (n - M) / numCoeffs;
						for (k = 1; k <= M; k++){
							value += 2.0 * freqResponce[k] * std::cos(x*k);
						}
						output[n] = value / numCoeffs;
					}
				} else {
					for (n = 0; n<numCoeffs; n++) {
						value = freqResponce[0];
						x = kTwoPi * (n - M) / numCoeffs;
						for (k = 1; k <= (numCoeffs / 2 - 1); k++){
							value += 2.0 * freqResponce[k] * std::cos(x*k);
						}
						output[n] = value / numCoeffs;
					}
				}
			} else {
				if (numCoeffs % 2) {
					for (n = 0; n<numCoeffs; n++) {
						value = 0;
						x = kTwoPi * (n - M) / numCoeffs;
						for (k = 1; k <= M; k++){
							value += 2.0 * freqResponce[k] * std::sin(x*k);
						}
						output[n] = value / numCoeffs;
					}
				} else {
					for (n = 0; n<numCoeffs; n++) {
						value = freqResponce[numCoeffs / 2] * std::sin(M_PI * (n - M));
						x = kTwoPi * (n - M) / numCoeffs;
						for (k = 1; k <= (numCoeffs / 2 - 1); k++){
							value += 2.0 * freqResponce[k] * std::sin(x*k);
						}
						output[n] = value / numCoeffs;
					}
				}
			}
		}
	};

	class DSP_SimpleLPF {
	public:
		struct Parameters {
			float g = 0.5f;
		};

		DSP_SimpleLPF() {}
		DSP_SimpleLPF(Parameters &parameters) : _params(parameters ) {}
		~DSP_SimpleLPF() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;
		}

		void Reset(int sampleRate) {
			_z1 = 0.0f;
		}

		float GetSample(float xn) {
			float yn  = _z1 * _params.g + xn * (1.0f - _params.g);
			_z1 = yn;

			return yn;
		}

	private:
		Parameters _params;
		float _z1 = 0.0f;
	};

	class DSP_SimpleAPF {
	public:
		enum FilterOrder {First, Second};
		struct Parameters {
			float frequency = 1000.0f;
			float q = 0.707f;
			float filterOrder = FilterOrder::First;
		};
		struct Coefficients {
			float a0 = 1.0f;
			float a1 = 0.0f;
			float a2 = 0.0f;
			float b1 = 0.0f;
			float b2 = 0.0f;
		};

		DSP_SimpleAPF() {}
		DSP_SimpleAPF(Parameters &parameters) : _params(parameters ) {}
		~DSP_SimpleAPF() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			if (_params.filterOrder == FilterOrder::First) {
				CalculateFirstOrderCoefficients();
			} else {
				CalculateSecondOrderCoefficients();
			}
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;

			_z1 = 0.0f;
			_z2 = 0.0f;

			SetParameters(_params);
		}

		float GetSample(float xn) {
			float yn = _coeffs.a0 * xn + _z1;

			_z1 = _coeffs.a1 * xn - _coeffs.b1 * yn + _z2;
			_z2 = _coeffs.a2 * xn - _coeffs.b2 * yn;

			return yn;
		}

	private:
		Parameters _params;
		Coefficients _coeffs;
		int _sampleRate = -1;
		float _z1 = 0.0f;
		float _z2 = 0.0f;

		void CalculateFirstOrderCoefficients() {
			float tanw0 = tan(M_PI * _params.frequency / _sampleRate);
			float alpha = (tanw0 - 1) / (tanw0 +1);

			_coeffs.a0 = alpha;
			_coeffs.a1 = 1.0f;
			_coeffs.a2 = 0.0f;
			_coeffs.b1 = alpha;
			_coeffs.b2 = 0.0f;
		}

		void CalculateSecondOrderCoefficients() {
			float bw = _params.frequency / _params.q;
			float tanw0 = tan(M_PI * bw / _sampleRate);
			float alpha = (tanw0 - 1) / (tanw0 +1);
			float beta = -cos(2.0f * M_PI * _params.frequency / _sampleRate);

			_coeffs.a0 = -alpha;
			_coeffs.a1 = beta * (1.0f - alpha);
			_coeffs.a2 = 1.0f;
			_coeffs.b1 = beta * (1.0f - alpha);
			_coeffs.b2 = -alpha;
		}
	};

	class DSP_SimpleDelay {
	public:
		enum DelayTimeUnit {MSec, Samples};
		struct Parameters {
			float delayTime = 0.0f;
			DelayTimeUnit delayTimeUnit = DelayTimeUnit::MSec;
			bool interpolate = false;
		};

		DSP_SimpleDelay() {}
		DSP_SimpleDelay(Parameters &parameters) : _params(parameters ) {}
		~DSP_SimpleDelay() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;

			if (_params.delayTimeUnit == DelayTimeUnit::MSec) {
				_delayInSamples = _params.delayTime * _sampleRate * 0.001;
			} else {
				_delayInSamples = _params.delayTime;
			}

			CreateDelayBuffer(1000);
		}

		void CreateDelayBuffer(float bufferLengthInMSec) {
			_bufferLength = bufferLengthInMSec * _sampleRate * 0.001;
			_buffer.CreateBuffer(_bufferLength);
		}

		float GetSample(float xn) {
			_buffer.WriteBuffer(xn);

			float yn = _buffer.ReadBuffer(_delayInSamples, _params.interpolate);

			return yn;
		}

		float ReadDelayAtSample(int sample) {
			return _buffer.ReadBuffer(sample);
		}

		float ReadDelayAtSample(float sample) {
			return _buffer.ReadBuffer(sample, _params.interpolate);
		}

		float ReadDelayAtPercentage(float percent) {
			return _buffer.ReadBuffer(percent * 0.01f * _bufferLength, _params.interpolate);
		}

		float ReadDelayAtTime(float mSec) {
			return _buffer.ReadBuffer(mSec * _sampleRate * 0.001f, _params.interpolate);
		}

	private:
		Parameters _params;
		CircularBuffer<float> _buffer;
		float _bufferLength;
		int _sampleRate;
		float _delayInSamples;
	};

	// Add RT60 calculation (new object?)
	class DSP_CombFilter {
	public:
		enum DelayTimeUnit {MSec, Samples};
		struct Parameters {
			float delayTime = 1.0f;
			DelayTimeUnit delayTimeUnit = DelayTimeUnit::MSec;
			bool interpolate = false;
			float g = 0.7f;
			float lpf_g = 0.0f;
		};

		DSP_CombFilter() {}
		DSP_CombFilter(Parameters &parameters) : _params(parameters ) {}
		~DSP_CombFilter() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			if (_params.delayTimeUnit == DelayTimeUnit::MSec) {
				_delayInSamples = _params.delayTime * _sampleRate * 0.001;
			} else {
				_delayInSamples = _params.delayTime;
			}
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;

			_z1 = 0.0f;

			SetParameters(_params);

			_buffer.CreateBuffer(_sampleRate);
		}

		float GetSample(float xn) {
			float yn = _buffer.ReadBuffer(_delayInSamples, _params.interpolate);
			float dn = (_z1 * _params.lpf_g) + yn;

			_z1 = dn;
			_buffer.WriteBuffer((dn * -_params.g) + xn);

			return yn;
		}

	private:
		Parameters _params;
		CircularBuffer<float> _buffer;
		int _sampleRate;
		float _delayInSamples;
		float _z1 = 0.0f;
	};

	class DSP_DelayAPF {
	public:
		enum DelayTimeUnit {MSec, Samples};
		struct Parameters {
			float delayTime = 1.0f;
			DelayTimeUnit delayTimeUnit = DelayTimeUnit::MSec;
			bool interpolate = true;
			float g = 0.7f;
		};

		DSP_DelayAPF() {}
		DSP_DelayAPF(Parameters &parameters) : _params(parameters ) {}
		~DSP_DelayAPF() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			if (_params.delayTimeUnit == DelayTimeUnit::MSec) {
				_delayInSamples = _params.delayTime * _sampleRate * 0.001;
			} else {
				_delayInSamples = _params.delayTime;
			}
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;

			SetParameters(_params);
			CreateDelayBuffer(1000);
		}

		float GetSample(float xn) {
			float yn = _params.g * xn + _buffer.ReadBuffer(_delayInSamples, _params.interpolate);

			_buffer.WriteBuffer(xn - _params.g * yn);

			return yn;
		}

	private:
		Parameters _params;
		int _sampleRate = -1;
		CircularBuffer<float> _buffer;
		float _delayInSamples;

		void CreateDelayBuffer(float bufferLengthInMSec) {
			int bufferLength = bufferLengthInMSec * (float)_sampleRate * 0.001f;
			_buffer.CreateBuffer(bufferLength);
		}
	};

	class DSP_SchroederReverb {
	public:
		struct Parameters {
			float cf1_delayTime = 45.0f;
			float cf1_g = 0.50f;
			float cf2_delayTime = 41.5f;
			float cf2_g = 0.62f;
			float cf3_delayTime = 37.2f;
			float cf3_g = 0.74f;
			float cf4_delayTime = 30.0f;
			float cf4_g = 0.86f;
			float apf1_delayTime = 5.0f;
			float apf1_g = 0.6f;
			float apf2_delayTime = 3.2f;
			float apf2_g = 0.5f;
			float apf3_delayTime = 2.75f;
			float apf3_g = 0.7f;
			float apf4_delayTime = 1.0f;
			float apf4_g = 0.65f;
			float dry = 0.3f;
		};

		DSP_SchroederReverb() {}
		DSP_SchroederReverb(Parameters &parameters) : _params(parameters ) {}
		~DSP_SchroederReverb() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			DSP_CombFilter::Parameters cfParams;

			cfParams.delayTime = _params.cf1_delayTime;
			cfParams.g = _params.cf1_g;
			_cf1.SetParameters(cfParams);

			cfParams.delayTime = _params.cf2_delayTime;
			cfParams.g = _params.cf2_g;
			_cf2.SetParameters(cfParams);

			cfParams.delayTime = _params.cf3_delayTime;
			cfParams.g = _params.cf3_g;
			_cf3.SetParameters(cfParams);

			cfParams.delayTime = _params.cf4_delayTime;
			cfParams.g = _params.cf4_g;
			_cf4.SetParameters(cfParams);

			DSP_DelayAPF::Parameters apfParams;

			apfParams.delayTime = _params.apf1_delayTime;
			apfParams.g = _params.apf1_g;
			_apf1.SetParameters(apfParams);

			apfParams.delayTime = _params.apf2_delayTime;
			apfParams.g = _params.apf2_g;
			_apf2.SetParameters(apfParams);

			apfParams.delayTime = _params.apf3_delayTime;
			apfParams.g = _params.apf3_g;
			_apf3.SetParameters(apfParams);

			apfParams.delayTime = _params.apf4_delayTime;
			apfParams.g = _params.apf4_g;
			_apf4.SetParameters(apfParams);
		}

		void Reset(int sampleRate) {
			_cf1.Reset(sampleRate);
			_cf2.Reset(sampleRate);
			_cf3.Reset(sampleRate);
			_cf4.Reset(sampleRate);
			_apf1.Reset(sampleRate);
			_apf2.Reset(sampleRate);
			_apf3.Reset(sampleRate);
			_apf4.Reset(sampleRate);
		}

		float GetSample(float xn) {
			float apf1 = _apf1.GetSample(xn);
			float apf2 = _apf2.GetSample(apf1);

			float cf = _cf1.GetSample(apf2)
			         + _cf2.GetSample(-apf2)
			         + _cf3.GetSample(apf2)
			         + _cf4.GetSample(-apf2);

			float apf3 = _apf3.GetSample(cf);
			float apf4 = _apf4.GetSample(apf3);

			float yn = xn * _params.dry + apf4 * (1.0f - _params.dry);

			return yn;
		}

	private:
		Parameters _params;
		DSP_CombFilter _cf1;
		DSP_CombFilter _cf2;
		DSP_CombFilter _cf3;
		DSP_CombFilter _cf4;
		DSP_DelayAPF _apf1;
		DSP_DelayAPF _apf2;
		DSP_DelayAPF _apf3;
		DSP_DelayAPF _apf4;
	};

	class DSP_LRCrossoverFilter {
	public:
		enum OutputMode {LowPass, HighPass};
		struct Parameters {
			float crossoverFrequency = 1000.0f;
			OutputMode outputMode = OutputMode::LowPass;
		};
		struct Output {
			float lp;
			float hp;
		};

		DSP_LRCrossoverFilter() {}
		DSP_LRCrossoverFilter(Parameters &parameters) : _params(parameters ) {}
		~DSP_LRCrossoverFilter() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			DSP_Biquad::Parameters filterParams;
			filterParams.frequency = _params.crossoverFrequency;
			
			filterParams.filterType = DSP_Biquad::FilterType::LinkwitzRileyLowPass;
			_lpf.SetParameters(filterParams);
			
			filterParams.filterType = DSP_Biquad::FilterType::LinkwitzRileyHighPass;
			_hpf.SetParameters(filterParams);
		}

		Output GetOutput() {return _output;}

		void Reset(int sampleRate) {
			_lpf.Reset(sampleRate);
			_hpf.Reset(sampleRate);

			SetParameters(_params);
		}

		float GetSample(float xn) {
			_output.lp = _lpf.GetSample(xn);
			_output.hp = -_hpf.GetSample(xn);

			float yn;
			switch (_params.outputMode) {
				case OutputMode::LowPass:
					yn = _output.lp;
					break;
				case OutputMode::HighPass:
					yn = _output.hp;
					break;
			}

			return yn;
		}

	private:
		Parameters _params;
		DSP_Biquad _lpf;
		DSP_Biquad _hpf;
		Output _output;
	};

	// Expander is broken
	class DSP_DynamicProcessor {
	public:
		enum DynamicProcessorType {Compressor, Limiter, Expander, Gate};
		struct Parameters {
			DynamicProcessorType dynamicProcessorType = DynamicProcessorType::Compressor;
			float threshold_dB = -12.0f;
			float ratio = 6.0f;
			float attackTime_mSec = 15.0f;
			float releaseTimes_mSec = 150.0f;
			float kneeWidth_dB = 3.0f;
			float makeupGain_dB = 3.0f;
		};

		DSP_DynamicProcessor() {}
		DSP_DynamicProcessor(Parameters &parameters) : _params(parameters ) {}
		~DSP_DynamicProcessor() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;

			DSP_AudioDetector::Parameters detectorParams;
			detectorParams.clampToUnityMax = false;
			detectorParams.detectDB = true;
			detectorParams.attackTime_mSec = _params.attackTime_mSec;
			detectorParams.releaseTime_mSec = _params.releaseTimes_mSec;
			_detector.SetParameters(detectorParams);
		}

		void Reset(int sampleRate) {
			_detector.Reset(sampleRate);

			SetParameters(_params);
		}

		float GetSample(float xn) {
			float detect_dB = _detector.GetSample(xn);

			float gainReduction = ComputeGain(detect_dB);

			float makeupGain = std::pow(10, _params.makeupGain_dB / 20.0f);

			return xn * gainReduction * makeupGain;
		}

	private:
		Parameters _params;
		DSP_AudioDetector _detector;

		float ComputeGain(float detect_dB) {
			float output_dB = 0.0f;

			switch(_params.dynamicProcessorType) {
				case DynamicProcessorType::Compressor:
					if (2.0f * (detect_dB - _params.threshold_dB) < -_params.kneeWidth_dB) { // Below threshold
						output_dB = detect_dB;
					} else if (2.0f * (std::abs(detect_dB - _params.threshold_dB)) <= _params.kneeWidth_dB) {
						output_dB = detect_dB + (((1.0f / _params.ratio) - 1.0f) * std::pow((detect_dB - _params.threshold_dB + (_params.kneeWidth_dB / 2.0f)), 2.0f)) / (2.0f * _params.kneeWidth_dB);
					} else if (2.0f * (detect_dB - _params.threshold_dB) > _params.kneeWidth_dB) {
						output_dB = _params.threshold_dB + (detect_dB - _params.threshold_dB) / _params.ratio;
					}
					break;
				case DynamicProcessorType::Limiter:
					if (2.0f * (detect_dB - _params.threshold_dB) < -_params.kneeWidth_dB) { // Below threshold
						output_dB = detect_dB;
					} else if (2.0f * (std::abs(detect_dB - _params.threshold_dB)) <= _params.kneeWidth_dB) {
						output_dB = detect_dB - std::pow((detect_dB - _params.threshold_dB + (_params.kneeWidth_dB / 2.0f)), 2.0f) / (2.0f * _params.kneeWidth_dB);
					} else if (2.0f * (detect_dB - _params.threshold_dB) > _params.kneeWidth_dB) {
						output_dB = _params.threshold_dB;
					}
					break;
				case DynamicProcessorType::Expander:
					if (2.0f * (detect_dB - _params.threshold_dB) > _params.kneeWidth_dB) { // Above threshold
						output_dB = detect_dB;
					} else if (2.0f * (std::abs(detect_dB - _params.threshold_dB)) > -_params.kneeWidth_dB) {
						output_dB = ((_params.ratio - 1.0f) * std::pow((detect_dB - _params.threshold_dB - (_params.kneeWidth_dB / 2.0f)), 2.0f)) / (2.0f * _params.kneeWidth_dB);
					} else if (2.0f * (detect_dB - _params.threshold_dB) <= -_params.kneeWidth_dB) {
						output_dB = _params.threshold_dB + (detect_dB - _params.threshold_dB) * _params.ratio;
					}
					break;
				case DynamicProcessorType::Gate:
					if (detect_dB >= _params.threshold_dB) { // Above threshold
						output_dB = detect_dB;
					} else {
						output_dB = -1.0e34;
					}
					break;
			}

			return std::pow(10, output_dB / 20.0f);
		}
	};

	class DSP_3BandCompressor {
	public:
		struct Parameters {};

		DSP_3BandCompressor() {}
		DSP_3BandCompressor(Parameters &parameters) : _params(parameters ) {}
		~DSP_3BandCompressor() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			DSP_LRCrossoverFilter::Parameters filterParams;
			filterParams.crossoverFrequency = 500.0f;
			_filter1.SetParameters(filterParams);
			filterParams.crossoverFrequency = 5000.0f;
			_filter2.SetParameters(filterParams);

			_params = parameters;
		}

		void Reset(int sampleRate) {
			_compressorL.Reset(sampleRate);
			_compressorM.Reset(sampleRate);
			_compressorH.Reset(sampleRate);
			_filter1.Reset(sampleRate);
			_filter2.Reset(sampleRate);

			SetParameters(_params);
		}

		float GetSample(float xn) {
			float sl = _filter2.GetSample(_filter1.GetSample(xn));
			float sm = _filter2.GetOutput().hp;
			float sh = -_filter2.GetOutput().hp;

			float yn = _compressorL.GetSample(sl)
			         + _compressorM.GetSample(sm)
			         + _compressorH.GetSample(sh);

			return yn;
		}

	private:
		Parameters _params;
		DSP_DynamicProcessor _compressorL;
		DSP_DynamicProcessor _compressorM;
		DSP_DynamicProcessor _compressorH;
		DSP_LRCrossoverFilter _filter1;
		DSP_LRCrossoverFilter _filter2;
	};
}