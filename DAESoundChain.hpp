#pragma once

#include "SoundChain.hpp"
#include "DEADSP.hpp"
#include <cmath>

namespace DEA {

	class BiquadSoundChain : public SoundChainBase {
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
			float frequency = 5000.0f;
			float q = 0.707f;
			float filterGain = 6.0f;
			FilterType filterType = FilterType::LowPass;
		};

		BiquadSoundChain() {}
		BiquadSoundChain(Parameters &parameters) : _params(parameters ) {}
		~BiquadSoundChain() {
			if (_x_z1) { delete[] _x_z1; }
			if (_x_z2) { delete[] _x_z2; }
			if (_y_z1) { delete[] _y_z1; }
			if (_y_z2) { delete[] _y_z2; }
		}

		Parameters GetParameters() {return _params;}
		void SetParameters(Parameters &parameters) {
			_params = parameters;
			CalculateCoefficients();
		}

	private:
		Parameters _params;
		float _a0 = 1.0f;
		float _a1 = 0.0f;
		float _a2 = 0.0f;
		float _b1 = 0.0f;
		float _b2 = 0.0f;
		float _c0 = 1.0f;
		float _d0 = 0.0f;
		float* _x_z1 = nullptr;
		float* _x_z2 = nullptr;
		float* _y_z1 = nullptr;
		float* _y_z2 = nullptr;

		void Reset() override {
			if (_x_z1) { delete[] _x_z1; }
			if (_x_z2) { delete[] _x_z2; }
			if (_y_z1) { delete[] _y_z1; }
			if (_y_z2) { delete[] _y_z2; }

			int numberOfChannels = ReadSettings().Channels;

			_x_z1 = new float[numberOfChannels];
			_x_z2 = new float[numberOfChannels];
			_y_z1 = new float[numberOfChannels];
			_y_z2 = new float[numberOfChannels];

			for (int channel = 0; channel < numberOfChannels; channel++) {
				_x_z1[channel] = 0.0f;
				_x_z2[channel] = 0.0f;
				_y_z1[channel] = 0.0f;
				_y_z2[channel] = 0.0f;
			}

			CalculateCoefficients();
		}

		void Process(float* buffPtr, int numberOfFrames) override {
			int numberOfChannels = ReadSettings().Channels;
			int numberOfSamples = numberOfFrames * numberOfChannels;

			for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
				for (int channel = 0; channel < numberOfChannels; channel++) {
					float xn = buffPtr[sample+channel];
					float yn = _a0 * xn
					         + _a1 * _x_z1[channel] 
					         + _a2 * _x_z2[channel] 
					         - _b1 * _y_z1[channel]
					         - _b2 * _y_z2[channel];

					_x_z2[channel] = _x_z1[channel];
					_x_z1[channel] = xn;

					_y_z2[channel] = _y_z1[channel];
					_y_z1[channel] = yn;

					buffPtr[sample+channel] = yn * _c0 + xn * _d0;
				}
			}
		}

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
			float w0 = 2.0f * M_PI * _params.frequency / ReadSettings().SampleRate;
			float cosw0 = cos(w0);
			float sinw0 = sin(w0);
			float d = 1.0f / _params.q;
			float beta = 0.5f * ((1.0f - d / 2.0f * sinw0) / (1.0f + d / 2.0f * sinw0));
			float gamma = (0.5f + beta) * cosw0;

			_a0 = (0.5f + beta - gamma) / 2.0f;
			_a1 = 0.5f + beta - gamma;
			_a2 = (0.5f + beta - gamma) / 2.0f;
			_b1 = -2.0f * gamma;
			_b2 = 2.0f * beta;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateHighPassCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / ReadSettings().SampleRate;
			float cosw0 = cos(w0);
			float sinw0 = sin(w0);
			float d = 1.0f / _params.q;
			float beta = 0.5f * ((1.0f - d / 2.0f * sinw0) / (1.0f + d / 2.0f * sinw0));
			float gamma = (0.5f + beta) * cosw0;

			_a0 = (0.5f + beta + gamma) / 2.0f;
			_a1 = -(0.5f + beta + gamma);
			_a2 = (0.5f + beta + gamma) / 2.0f;
			_b1 = -2.0f * gamma;
			_b2 = 2.0f * beta;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateBandPassCoefficients() {
			float k = tan((M_PI * _params.frequency) / ReadSettings().SampleRate);
			float k2 = k*k;
			float delta = k2 * _params.q + k + _params.q;

			_a0 = k / _params.q;
			_a1 = 0.0f;
			_a2 = -k / delta;
			_b1 = (2.0f * _params.q * (k2 - 1)) / delta;
			_b2 = (k2 * _params.q - k + _params.q) / delta;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateBandStopCoefficients() {
			float k = tan((M_PI * _params.frequency) / ReadSettings().SampleRate);
			float k2 = k*k;
			float delta = k2 * _params.q + k + _params.q;

			_a0 = (_params.q * (k2 + 1)) / delta;
			_a1 = (2.0f * _params.q * (k2 - 1)) / delta;
			_a2 = (_params.q * (k2 + 1)) / delta;
			_b1 = (2.0f * _params.q * (k2 - 1)) / delta;
			_b2 = (k2 * _params.q - k + _params.q) / delta;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateLowShelfCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / ReadSettings().SampleRate;
			float mu = std::pow(10, _params.filterGain / 20.0f);
			float beta = 4.0f / (1.0f + mu);
			float delta = beta * tan(w0 / 2.0f);
			float gamma = (1.0f - delta) / (1.0f + delta);

			_a0 = (1 - gamma) / 2.0f;
			_a1 = (1 - gamma) / 2.0f;
			_a2 = 0.0f;
			_b1 = -gamma;
			_b2 = 0.0f;
			_c0 = mu - 1.0f;
			_d0 = 1.0f;
		}

		void CalculateHighShelfCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / ReadSettings().SampleRate;
			float mu = std::pow(10, _params.filterGain / 20.0f);
			float beta = (1.0f + mu) / 4.0f;
			float delta = beta * tan(w0 / 2.0f);
			float gamma = (1.0f - delta) / (1.0f + delta);

			_a0 = (1 + gamma) / 2.0f;
			_a1 = (1 + gamma) / -2.0f;
			_a2 = 0.0f;
			_b1 = -gamma;
			_b2 = 0.0f;
			_c0 = mu - 1.0f;
			_d0 = 1.0f;
		}

		void CalculateButterworthLowPassCoefficients() {
			float c = 1.0f / tan((M_PI * _params.frequency) / ReadSettings().SampleRate);

			_a0 = 1.0f / (1.0f + std::sqrt(2.0f * c) + c*c);
			_a1 = 2.0f * _a0;
			_a2 = _a0;
			_b1 = 2.0f * _a0 * (1.0f - c*c);
			_b2 = _a0 * (1.0f - std::sqrt(2.0f * c) + c*c);
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateButterworthHighPassCoefficients() {
			float c = tan((M_PI * _params.frequency) / ReadSettings().SampleRate);

			_a0 = 1.0f / (1.0f + std::sqrt(2.0f * c) + c*c);
			_a1 = -2.0f * _a0;
			_a2 = _a0;
			_b1 = 2.0f * _a0 * (c*c - 1.0f);
			_b2 = _a0 * (1.0f - std::sqrt(2.0f * c) + c*c);
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateLinkwitzRileyLowPassCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / ReadSettings().SampleRate;
			float omega = M_PI * _params.frequency;
			float kappa = omega / tan(w0);
			float delta = kappa*kappa + omega*omega + 2.0f * kappa * omega;
			float omega2 = omega*omega;
			float kappa2 = kappa*kappa;

			_a0 = omega2 / delta;
			_a1 = 2.0f * (omega2 / delta);
			_a2 = omega2 / delta;
			_b1 = (-2.0f * kappa2 + 2.0f * omega2) / delta;
			_b2 = (-2.0f * kappa * omega + kappa2 + omega2) / delta;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateLinkwitzRileyHighPassCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / ReadSettings().SampleRate;
			float omega = M_PI * _params.frequency;
			float kappa = omega / tan(w0);
			float delta = kappa*kappa + omega*omega + 2.0f * kappa * omega;
			float omega2 = omega*omega;
			float kappa2 = kappa*kappa;

			_a0 = kappa2 / delta;
			_a1 = (-2.0f * kappa2) / delta;
			_a2 = kappa2 / delta;
			_b1 = (-2.0f * kappa2 + 2.0f * omega2) / delta;
			_b2 = (-2.0f * kappa * omega + kappa2 + omega2) / delta;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateAllPass1stOrderCoefficients() {
			float w0 = M_PI * _params.frequency / ReadSettings().SampleRate;
			float tanw0 = tan(w0);
			float alpha = (tanw0 - 1) / (tanw0 +1);

			_a0 = alpha;
			_a1 = 1.0f;
			_a2 = 0.0f;
			_b1 = alpha;
			_b2 = 0.0f;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateAllPass2ndOrderCoefficients() {
			float bw = _params.frequency / _params.q;
			float tanbwpisr = tan(bw * M_PI / ReadSettings().SampleRate);
			float alpha = (tanbwpisr - 1) / (tanbwpisr +1);
			float beta = -cos(2.0f * M_PI * _params.frequency / ReadSettings().SampleRate);

			_a0 = -alpha;
			_a1 = beta * (1.0f - alpha);
			_a2 = 1.0f;
			_b1 = beta * (1.0f - alpha);
			_b2 = -alpha;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateLowPass1stOrderCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / ReadSettings().SampleRate;
			float cosw0 = cos(w0);
			float sinw0 = sin(w0);
			float gamma = cosw0 / (1.0f + sinw0);

			_a0 = (1.0f - gamma) / 2.0f;
			_a1 = (1.0f - gamma) / 2.0f;
			_a2 = 0.0f;
			_b1 = -gamma;
			_b2 = 0.0f;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculateHighPass1stOrderCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / ReadSettings().SampleRate;
			float cosw0 = cos(w0);
			float sinw0 = sin(w0);
			float gamma = cosw0 / (1.0f + sinw0);

			_a0 = (1.0f + gamma) / 2.0f;
			_a1 = -((1.0f + gamma) / 2.0f);
			_a2 = 0.0f;
			_b1 = -gamma;
			_b2 = 0.0f;
			_c0 = 1.0f;
			_d0 = 0.0f;
		}

		void CalculatePeakCoefficients() {
			float w0 = 2.0f * M_PI * _params.frequency / ReadSettings().SampleRate;
			float mu = std::pow(10, _params.filterGain / 20.0f);
			float zeta = 4.0f / (1.0f + mu);
			float zetatanw02q = zeta * tan(w0 / 2.0f * _params.q);
			float beta = 0.5f * (1.0f - zetatanw02q) / (1.0f + zetatanw02q);
			float gamma = (0.5f + beta) * cos(w0);

			_a0 = 0.5f - beta;
			_a1 = 0.0f;
			_a2 = -(0.5f - beta);
			_b1 = -2.0f * gamma;
			_b2 =2.0f * beta;
			_c0 = mu - 1.0f;
			_d0 = 1.0f;
		}
	};

	class SimpleResonatorSoundChain : public SoundChainBase {
	public:
		SimpleResonatorSoundChain(float frequency) : _f0(frequency) {}
		~SimpleResonatorSoundChain() {
			if (_y_z1) { delete[] _y_z1; }
			if (_y_z2) { delete[] _y_z2; }
		}

		void SetFrequency(float frequency) {
			_f0 = frequency;
			CalculateCoefficients();
		}
		void SetQ(float q) {
			_q = q;
			CalculateCoefficients();
		}
		void SetFilterSettings(float frequency, float q) {
			_f0 = frequency;
			_q = q;
			CalculateCoefficients();
		}

	private:
		float _a0 = 1.0f;
		float _b1 = 0.0f;
		float _b2 = 0.0f;
		float* _y_z1 = nullptr;
		float* _y_z2 = nullptr;
		int _fs = 0;
		float _f0 = 5000.0f;
		float _q = 1.0f;

		void Reset() override {
			if (_y_z1) { delete[] _y_z1; }
			if (_y_z2) { delete[] _y_z2; }

			_y_z1 = new float[ReadSettings().Channels];
			_y_z2 = new float[ReadSettings().Channels];

			int numberOfChannels = ReadSettings().Channels;
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_y_z1[channel] = 0.0f;
				_y_z2[channel] = 0.0f;
			}

			_fs = ReadSettings().SampleRate;
			CalculateCoefficients();
		}

		void Process(float* buffPtr, int numberOfFrames) override {
			int numberOfChannels = ReadSettings().Channels;
			int numberOfSamples = numberOfFrames * numberOfChannels;

			for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
				for (int channel = 0; channel < numberOfChannels; channel++) {
					float xn = buffPtr[sample+channel];
					float yn = _a0 * xn - _b1 * _y_z1[channel] - _b2 * _y_z2[channel];

					_y_z2[channel] = _y_z1[channel];
					_y_z1[channel] = yn;

					buffPtr[sample+channel] = yn;
				}
			}
		}

		void CalculateCoefficients() {
			float w0 = 2.0f * M_PI * _f0 / _fs;
			float cosw0 = cos(w0);
			float bw = _f0 / _q;

			_b2 = std::exp(-2.0f * M_PI * (bw / _fs));
			_b1 = ((-4.0f * _b2) / (1.0f + _b2)) * cosw0;
			_a0 = (1.0f - _b2) * std::sqrt(1.0f - _b1*_b1 / 4.0f * _b2);
		}
	};

	class SmithAngellResonatorSoundChain : public SoundChainBase {
	public:
		SmithAngellResonatorSoundChain(float frequency) : _f0(frequency) {}
		~SmithAngellResonatorSoundChain() {
			if (_x_z1) { delete[] _x_z1; }
			if (_x_z2) { delete[] _x_z2; }
			if (_y_z1) { delete[] _y_z1; }
			if (_y_z2) { delete[] _y_z2; }
		}

		void SetFrequency(float frequency) {
			_f0 = frequency;
			CalculateCoefficients();
		}
		void SetQ(float q) {
			_q = q;
			CalculateCoefficients();
		}
		void SetFilterSettings(float frequency, float q) {
			_f0 = frequency;
			_q = q;
			CalculateCoefficients();
		}

	private:
		float _a0 = 1.0f;
		float _a2 = 0.0f;
		float _b1 = 0.0f;
		float _b2 = 0.0f;
		float* _x_z1 = nullptr;
		float* _x_z2 = nullptr;
		float* _y_z1 = nullptr;
		float* _y_z2 = nullptr;
		int _fs = 0;
		float _f0 = 5000.0f;
		float _q = 1.0f;

		void Reset() override {
			if (_x_z1) { delete[] _x_z1; }
			if (_x_z2) { delete[] _x_z2; }
			if (_y_z1) { delete[] _y_z1; }
			if (_y_z2) { delete[] _y_z2; }

			_x_z1 = new float[ReadSettings().Channels];
			_x_z2 = new float[ReadSettings().Channels];
			_y_z1 = new float[ReadSettings().Channels];
			_y_z2 = new float[ReadSettings().Channels];

			int numberOfChannels = ReadSettings().Channels;
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_x_z1[channel] = 0.0f;
				_x_z2[channel] = 0.0f;
				_y_z1[channel] = 0.0f;
				_y_z2[channel] = 0.0f;
			}

			_fs = ReadSettings().SampleRate;
			CalculateCoefficients();
		}

		void Process(float* buffPtr, int numberOfFrames) override {
			int numberOfChannels = ReadSettings().Channels;
			int numberOfSamples = numberOfFrames * numberOfChannels;

			for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
				for (int channel = 0; channel < numberOfChannels; channel++) {
					float xn = buffPtr[sample+channel];
					float yn = _a0 * xn 
					         + _a2 * _x_z2[channel] 
					         - _b1 * _y_z1[channel]
					         - _b2 * _y_z2[channel];

					_x_z2[channel] = _x_z1[channel];
					_x_z1[channel] = xn;

					_y_z2[channel] = _y_z1[channel];
					_y_z1[channel] = yn;

					buffPtr[sample+channel] = yn;
				}
			}
		}

		void CalculateCoefficients() {
			float w0 = 2.0f * M_PI * _f0 / _fs;
			float cosw0 = cos(w0);
			float bw = _f0 / _q;

			_b2 = std::exp(-2.0f * M_PI * (bw / _fs));
			_b1 = ((-4.0f * _b2) / (1.0f + _b2)) * cosw0;
			_a0 = 1.0f - std::sqrt(_b2);
			_a2 = -_a0;
		}
	};

	template <class DSP_T>
	class DSPSoundChain : public SoundChainBase {
	public:
		DSPSoundChain() {}
		DSPSoundChain(typename DSP_T::Parameters &parameters) : _params(parameters ) {}
		~DSPSoundChain() {
			if (_bank) { delete[] _bank; }
		}

		typename DSP_T::Parameters GetParameters() {return _params;}
		void SetParameters(typename DSP_T::Parameters &parameters) {
			_params = parameters;
			PropogateParameters();
		}

	private:
		typename DSP_T::Parameters _params;
		DSP_T* _bank = nullptr;

		void Reset() override {
			if (_bank) { delete[] _bank; }

			int numberOfChannels = ReadSettings().Channels;
			int sampleRate = ReadSettings().SampleRate;

			_bank = new DSP_T[numberOfChannels];
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].Reset(sampleRate);
			}

			PropogateParameters();
		}

		void Process(float* buffPtr, int numberOfFrames) override {
			int numberOfChannels = ReadSettings().Channels;
			int numberOfSamples = numberOfFrames * numberOfChannels;

			for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
				for (int channel = 0; channel < numberOfChannels; channel++) {
					buffPtr[sample+channel] = _bank[channel].GetSample(buffPtr[sample+channel]);
				}
			}
		}

		void PropogateParameters() {
			if (!_bank) return;

			int numberOfChannels = ReadSettings().Channels;
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].SetParameters(_params);
			}
		}
	};

	// TODO: pingpong not pingponging
	template <class DSP_T>
	class DSPPingPongSoundChain : public SoundChainBase {
	public:
		DSPPingPongSoundChain() {}
		DSPPingPongSoundChain(typename DSP_T::Parameters &parameters) : _params(parameters ) {}
		~DSPPingPongSoundChain() {
			if (_bank) { delete[] _bank; }
			if (_outSample) { delete[] _outSample; }
		}

		typename DSP_T::Parameters GetParameters() {return _params;}
		void SetParameters(typename DSP_T::Parameters &parameters) {
			_params = parameters;
			
			PropogateParameters();
		}

	private:
		typename DSP_T::Parameters _params;
		DSP_T* _bank = nullptr;
		float* _outSample = nullptr;

		void Reset() override {
			if (_bank) { delete[] _bank; }
			if (_outSample) { delete[] _outSample; }

			int numberOfChannels = ReadSettings().Channels;
			int sampleRate = ReadSettings().SampleRate;

			_bank = new DSP_T[numberOfChannels];
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].Reset(sampleRate);
			}
			_outSample = new float[numberOfChannels];
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_outSample[channel] = 0.0f;
			}

			PropogateParameters();
		}

		void Process(float* buffPtr, int numberOfFrames) override {
			int numberOfChannels = ReadSettings().Channels;
			int numberOfSamples = numberOfFrames * numberOfChannels;

			for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
				_outSample[0] = _bank[0].GetSample(buffPtr[sample+numberOfChannels-1]);
				for (int channel = 1; channel < numberOfChannels; channel++) {
					_outSample[channel] = _bank[channel].GetSample(buffPtr[sample+channel-1]);
				}
				for (int channel = 0; channel < numberOfChannels; channel++) {
					buffPtr[sample+channel] = _outSample[channel];
				}
			}
		}

		void PropogateParameters() {
			if (!_bank) return;

			int numberOfChannels = ReadSettings().Channels;
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].SetParameters(_params);
			}
		}
	};

	template <class DSP_T>
	class DSPLCRSoundChain : public SoundChainBase {
	public:
		DSPLCRSoundChain() {}
		DSPLCRSoundChain(typename DSP_T::Parameters &parameters) : _params(parameters ) {}
		~DSPLCRSoundChain() {
			if (_bank) { delete[] _bank; }
		}

		typename DSP_T::Parameters GetParameters() {return _params;}
		void SetParameters(typename DSP_T::Parameters &parameters) {
			_params = parameters;
			
			PropogateParameters();
		}

	private:
		typename DSP_T::Parameters _params;
		DSP_T* _bank = nullptr;

		void Reset() override {
			if (_bank) { delete[] _bank; }

			int numberOfChannels = ReadSettings().Channels;
			int sampleRate = ReadSettings().SampleRate;

			_bank = new DSP_T[numberOfChannels];
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].Reset(sampleRate);
			}

			PropogateParameters();
		}

		void Process(float* buffPtr, int numberOfFrames) override {
			int numberOfChannels = ReadSettings().Channels;
			int numberOfSamples = numberOfFrames * numberOfChannels;

			for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
				float sum = 0;
				for (int channel = 0; channel < numberOfChannels; channel++) {
					sum += buffPtr[sample+channel];
				}
				sum *= 0.5f;
				for (int channel = 0; channel < numberOfChannels; channel++) {
					buffPtr[sample+channel] = _bank[channel].GetSample(buffPtr[sample+channel] + sum);
				}
			}
		}

		void PropogateParameters() {
			if (!_bank) return;

			int numberOfChannels = ReadSettings().Channels;
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].SetParameters(_params);
			}
		}
	};

	class PingPongDelaySoundChain : public SoundChainBase {
	public:
		PingPongDelaySoundChain() {}
		PingPongDelaySoundChain(DSP_AnalogDelay::Parameters &parameters) : _params(parameters ) {}
		~PingPongDelaySoundChain() {
			if (_bank) { delete[] _bank; }
			if (_outSample) { delete[] _outSample; }
		}

		DSP_AnalogDelay::Parameters GetParameters() {return _params;}
		void SetParameters(DSP_AnalogDelay::Parameters &parameters) {
			_params = parameters;
			
			PropogateParameters();
		}

	private:
		DSP_AnalogDelay::Parameters _params;
		DSP_AnalogDelay* _bank = nullptr;
		float* _outSample = nullptr;

		void Reset() override {
			if (_bank) { delete[] _bank; }
			if (_outSample) { delete[] _outSample; }

			int numberOfChannels = ReadSettings().Channels;
			int sampleRate = ReadSettings().SampleRate;

			_bank = new DSP_AnalogDelay[numberOfChannels];
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].Reset(sampleRate);
			}
			_outSample = new float[numberOfChannels];
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_outSample[channel] = 0.0f;
			}

			PropogateParameters();
		}

		void Process(float* buffPtr, int numberOfFrames) override {
			int numberOfChannels = ReadSettings().Channels;
			int numberOfSamples = numberOfFrames * numberOfChannels;

			for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
				_outSample[0] = _bank[0].GetSample(buffPtr[sample+numberOfChannels-1]);
				for (int channel = 1; channel < numberOfChannels; channel++) {
					_outSample[channel] = _bank[channel].GetSample(buffPtr[sample+channel-1]);
				}
				for (int channel = 0; channel < numberOfChannels; channel++) {
					buffPtr[sample+channel] = _outSample[channel];
				}
			}
		}

		void PropogateParameters() {
			if (!_bank) return;

			int numberOfChannels = ReadSettings().Channels;
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].SetParameters(_params);
			}
		}
	};

	class LCRDelaySoundChain : public SoundChainBase {
	public:
		LCRDelaySoundChain() {}
		LCRDelaySoundChain(DSP_AnalogDelay::Parameters &parameters) : _params(parameters ) {}
		~LCRDelaySoundChain() {
			if (_bank) { delete[] _bank; }
		}

		DSP_AnalogDelay::Parameters GetParameters() {return _params;}
		void SetParameters(DSP_AnalogDelay::Parameters &parameters) {
			_params = parameters;
			
			PropogateParameters();
		}

	private:
		DSP_AnalogDelay::Parameters _params;
		DSP_AnalogDelay* _bank = nullptr;

		void Reset() override {
			if (_bank) { delete[] _bank; }

			int numberOfChannels = ReadSettings().Channels;
			int sampleRate = ReadSettings().SampleRate;

			_bank = new DSP_AnalogDelay[numberOfChannels];
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].Reset(sampleRate);
			}

			PropogateParameters();
		}

		void Process(float* buffPtr, int numberOfFrames) override {
			int numberOfChannels = ReadSettings().Channels;
			int numberOfSamples = numberOfFrames * numberOfChannels;

			for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
				float sum = 0;
				for (int channel = 0; channel < numberOfChannels; channel++) {
					sum += buffPtr[sample+channel];
				}
				sum *= 0.5f;
				for (int channel = 0; channel < numberOfChannels; channel++) {
					buffPtr[sample+channel] = _bank[channel].GetSample(buffPtr[sample+channel] + sum);
				}
			}
		}

		void PropogateParameters() {
			if (!_bank) return;

			int numberOfChannels = ReadSettings().Channels;
			for (int channel = 0; channel < numberOfChannels; channel++) {
				_bank[channel].SetParameters(_params);
			}
		}
	};
}