#pragma once

#include <memory>
#include <cmath>

namespace DEA {

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
			float filterGain = 6.0f;
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
		struct Parameters {
			float frequency = 1000.0f;
			float q = 0.707f;
			FilterOrder filterOrder = FilterOrder::Second;
		};
		struct Output {
			float lpf;
			float apf;
			float hpf;
		};

		DSP_MutiFilter() {
			SetParameters(_params);
		}
		DSP_MutiFilter(Parameters &parameters) : _params(parameters ) {
			SetParameters(_params);
		}
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

		void Reset(int sampleRate) {
			_lpf.Reset(sampleRate);
		}

		Output GetSample(float xn) {
			Output output;

			output.lpf = _lpf.GetSample(xn);
			output.hpf = xn - output.lpf;
			output.apf = output.lpf - output.hpf;

			return output;
		}

	private:
		Parameters _params;
		DSP_Biquad _lpf;
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
		enum AudioDetectMode { Peak, MS, RMS};
		struct Parameters {
			float attackTime = 0.0f;
			float releaseTime = 0.0f;
			AudioDetectMode detectMode = AudioDetectMode::Peak;
			bool detectDB = false;
			bool clampToUnityMax = true;
		};

		DSP_AudioDetector() {}
		DSP_AudioDetector(Parameters &parameters) : _params(parameters ) {}
		~DSP_AudioDetector() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			if (_params.attackTime != parameters.attackTime) {
				_attackTime = exp(ANAOG_TC / (parameters.attackTime * _sampleRate * 0.001f));
			}
			if (_params.releaseTime != parameters.releaseTime) {
				_releaseTime = exp(ANAOG_TC / (parameters.releaseTime * _sampleRate * 0.001f));
			}
			
			_params = parameters;
		}

		void Reset(int sampleRate) {
			_sampleRate = sampleRate;
			_lastValue = 0.0f;
		}

		float GetSample(float xn) {
			float input = xn;

			// Rectifier
			input = abs(input);
			if (_params.detectMode == AudioDetectMode::MS || _params.detectMode == AudioDetectMode::RMS) {
				input *= input;
			}

			// RC Simulation
			float currentValue;
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
				currentValue = pow(currentValue, 0.5f);
			}
			if (_params.detectDB) {
				currentValue = 20.0f * log10(currentValue);
			}

			return currentValue;
		}

	private:
		const float ANAOG_TC = log(36.7f);

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
			float attackTime = 10.0f;
			float releaseTime = 10.0f;
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
			detectorParams.attackTime = _params.attackTime;
			detectorParams.releaseTime = _params.releaseTime;
			detectorParams.detectMode = DSP_AudioDetector::AudioDetectMode::RMS;
			detectorParams.detectDB = true;
			detectorParams.clampToUnityMax = false;
			_detector.SetParameters(detectorParams);
		}

		void Reset(int sampleRate) {
			_filter.Reset(sampleRate);
			_detector.Reset(sampleRate);
			_maxFilterFrequency = sampleRate / 2.0f;
		}

		float GetSample(float xn) {
			float thresholdValue = pow(10.0f, _params.threshold / 20.0f);

			float detectDB = _detector.GetSample(xn);
			float detectValue = pow(10.0f, detectDB / 20.0f);
			float deltaValue = detectValue - thresholdValue;

			auto filterParams = _filter.GetParameters();
			filterParams.frequency = _params.frequency;

			if (deltaValue > 0.0f) {
				float modulatorValue = deltaValue * _params.sensitivity;
				filterParams.frequency = modulatorValue * (_maxFilterFrequency - filterParams.frequency) + filterParams.frequency;
			}

			_filter.SetParameters(filterParams);

			return _filter.GetSample(xn);
		}

	private:
		Parameters _params;
		DSP_Biquad _filter;
		DSP_AudioDetector _detector;
		float _maxFilterFrequency = -1;
	};

	class DSP_Phaser {
	public:
		struct Parameters {
			float lfoRate;
			float lfoDepth;
			float intensity;
		};

		DSP_Phaser() {}
		DSP_Phaser(Parameters &parameters) : _params(parameters ) {}
		~DSP_Phaser() {}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;
		}

		void Reset(int sampleRate) {
			for (int i = 0; i < STAGES; i++) {
				_filterBank[i].filter.Reset(sampleRate);

				float scale = (i / STAGES) * (i / STAGES);
				_filterBank[i].minimumFrequency = 200.0f * scale;
				_filterBank[i].maximumFrequency = 20000.0f * scale;
			}

			_lfo.Reset(sampleRate);
			auto lfoParams = _lfo.GetParameters();
			lfoParams.frequency = _params.lfoRate;
			lfoParams.waveformType = DSP_LFO::WaveformType::Sin;
			_lfo.SetParameters(lfoParams);
		}

		float GetSample(float xn) {
			auto lfo = _lfo.GetSample();
			float lfoValue = (lfo.normalOutput * 0.5f + 0.5f) * _params.lfoDepth;
			for (int i = 0; i < STAGES; i++) {
				_filterBank[i].SetFrequencyFromScale(lfoValue);
			}

			float gamma = _filterBank[0].filter.GetCoefficients().a0;
			for (int i = 1; i < STAGES; i++) {
				gamma *= _filterBank[i].filter.GetCoefficients().a0;
			}
			float alpha = 1.0f / (1.0f + _params.intensity * gamma);

			// Missing feedback loop, need the S value
			float sn = 0.0f;
			// for (int i = 0; i < STAGES; i++) {
			// }

			float u = alpha * (xn - _params.intensity * sn);

			float yn = _filterBank[0].filter.GetSample(u);
			for (int i = 1; i < STAGES; i++) {
				yn = _filterBank[i].filter.GetSample(yn);
			}

			return xn * 0.707f + yn * 0.707f ;
		}

	private:
		struct AllPassFilter {
			DSP_Biquad filter;
			DSP_Biquad::Parameters filterParams;
			float minimumFrequency = 16.0f;
			float maximumFrequency = 20000.0f;

			AllPassFilter() {
				filterParams = filter.GetParameters();
				filterParams.filterType = DSP_Biquad::FilterType::AllPass1stOrder;
				filter.SetParameters(filterParams);
			}

			void SetFrequencyFromScale(float scale) {
				float halfRange = (maximumFrequency - minimumFrequency) / 2.0f;
				float midPoint = halfRange + minimumFrequency;
				filterParams.frequency = scale * halfRange + midPoint;
				filter.SetParameters(filterParams);
			}
		};

		const int static STAGES = 6;
		Parameters _params;
		AllPassFilter _filterBank[STAGES];
		DSP_LFO _lfo;
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
			return _buffer[_writeIndex];
		}

		T ReadBuffer(float delayInFractionalSamples, bool interpolate = true) {
			T start = ReadBuffer((int)delayInFractionalSamples);
			if (!interpolate) return start;

			T end = ReadBuffer((int)delayInFractionalSamples + 1);
			float t = delayInFractionalSamples - (int) delayInFractionalSamples;

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
			_samplesPerMS = _sampleRate / 1000.0f;

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
		float _samplesPerMS;
	};
}