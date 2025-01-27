#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

struct SoundChainSettings {
	int SampleRate;
	int Channels;
};

class SoundChainBase {
protected:
	virtual void Reset() {};
	virtual void Process(float* buffPtr, int numberOfFrames) {};

public:
	void Initialize(SoundChainSettings &soundChainSettings) {
		_settings = soundChainSettings;
		
		if (_previous) {
			_previous->Initialize(soundChainSettings);
		}

		Reset();
		_initialized = true;
	}

	bool IsInitialized() {return _initialized;}
	SoundChainSettings ReadSettings() {return _settings;}

	SoundChainBase* GetPrevious() {return _previous;}
	void SetPrevious(SoundChainBase* previous) {_previous = previous;}

	void ReadSamples(float* buffPtr, int numberOfFrames) {
		if (!_initialized) return;

		if (_previous) {
			_previous->ReadSamples(buffPtr, numberOfFrames);
		} else {
			int numberOfSamples = numberOfFrames * _settings.Channels;
			for (int sample = 0; sample < numberOfSamples; sample++) {
				buffPtr[sample] = 0.0f;
			}
		}

		Process(buffPtr, numberOfFrames);
	}

private:
	SoundChainSettings _settings;
	SoundChainBase* _previous = nullptr;
	bool _initialized = false;
};

class TrimSoundChain : public SoundChainBase {
public:
	struct Parameters {
		float amplitude = 1.0f;
	};

	TrimSoundChain() {}
	TrimSoundChain(Parameters &parameters) : _params(parameters ) {}

	Parameters GetParameters() {return _params;}
	void SetParameters(Parameters &parameters) {
		_params = parameters;
	}

private:
	Parameters _params;

	void Process(float* buffPtr, int numberOfFrames) override {
		int numberOfSamples = numberOfFrames * ReadSettings().Channels;
		for (int sample = 0; sample < numberOfSamples; sample++) {
			buffPtr[sample] = _params.amplitude * buffPtr[sample];
		}
	}
};

class SinewaveSoundChain : public SoundChainBase {
public:
	struct Parameters {
		double frequency = 440.0;
		float amplitude = 0.25f;
	};

	SinewaveSoundChain() {}
	SinewaveSoundChain(Parameters &parameters) : _params(parameters ) {}
	SinewaveSoundChain(double frequency) {
		_params.frequency = frequency;
	}

	Parameters GetParameters() {return _params;}
	void SetParameters(Parameters &parameters) {
		_params = parameters;
		_stepSize = _params.frequency / ReadSettings().SampleRate;
	}

private:
	Parameters _params;
	double _phase = 0.0;
	double _stepSize = 0.0;

	void Reset() override {
		_stepSize = _params.frequency / ReadSettings().SampleRate;
	}

	void Process(float* buffPtr, int numberOfFrames) override {
		int numberOfChannels = ReadSettings().Channels;
		int numberOfSamples = numberOfFrames * numberOfChannels;
		for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
			float value =  _params.amplitude * sin(_phase * 2 * M_PI);

			for (int channel = 0; channel < numberOfChannels; channel++) {
				buffPtr[sample + channel] += value;
			}

			_phase += _stepSize;
			if (_phase > 1) _phase -= 1;
		}
	}
};

class WhiteNoiseSoundChain : public SoundChainBase {
public:
	struct Parameters {
		float amplitude = 0.2f;
	};

	WhiteNoiseSoundChain() {}
	WhiteNoiseSoundChain(Parameters &parameters) : _params(parameters ) {}

	Parameters GetParameters() {return _params;}
	void SetParameters(Parameters &parameters) {
		_params = parameters;
	}

private:
	Parameters _params;

	void Process(float* buffPtr, int numberOfFrames) {
		int numberOfChannels = ReadSettings().Channels;
		int numberOfSamples = numberOfFrames * numberOfChannels;
		for (int sample = 0; sample < numberOfSamples;  sample += numberOfChannels) {
			for (int channel = 0; channel < numberOfChannels; channel++) {
				buffPtr[sample + channel] += _params.amplitude * ((rand() + 0.0f)/RAND_MAX * 2 - 1);
			}
		}
	}
};

class CookEQSoundChain : public SoundChainBase {
public:

	enum FilterType {LowPass, HighPass, BandPass, Notch, AllPass};
	struct Parameters {
		float frequency = 5000.0f;
		float q = 0.707f;
		FilterType filterType = FilterType::LowPass;
	};

	CookEQSoundChain() {}
	CookEQSoundChain(Parameters &parameters) : _params(parameters ) {}
	~CookEQSoundChain() {
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
	float _b0 = 0.0f;
	float _b1 = 0.0f;
	float _b2 = 0.0f;
	float* _x_z1 = nullptr;
	float* _x_z2 = nullptr;
	float* _y_z1 = nullptr;
	float* _y_z2 = nullptr;
	int _fs = -1;
	FilterType _filterType = FilterType::LowPass;

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

		_fs = ReadSettings().SampleRate;
		CalculateCoefficients();
	}

	void Process(float* buffPtr, int numberOfFrames) override {
		int numberOfChannels = ReadSettings().Channels;
		int numberOfSamples = numberOfFrames * numberOfChannels;

		for (int sample = 0; sample < numberOfSamples; sample += numberOfChannels) {
			for (int channel = 0; channel < numberOfChannels; channel++) {
				float xn = buffPtr[sample+channel];
				float yn = (_b0/_a0) * xn 
				         + (_b1/_a0) * _x_z1[channel] 
				         + (_b2/_a0) * _x_z2[channel] 
				         - (_a1/_a0) * _y_z1[channel]
				         - (_a2/_a0) * _y_z2[channel];

				_x_z2[channel] = _x_z1[channel];
				_x_z1[channel] = xn;

				_y_z2[channel] = _y_z1[channel];
				_y_z1[channel] = yn;

				buffPtr[sample+channel] = yn;
			}
		}
	}

	void CalculateCoefficients() {
		switch (_params.filterType) {
			case FilterType::LowPass:
				LowPassCoefficients();
				break;
			case FilterType::HighPass:
				HighPassCoefficients();
				break;
			case FilterType::BandPass:
				BandPassCoefficients();
				break;
			case FilterType::Notch:
				NotchCoefficients();
				break;
			case FilterType::AllPass:
				AllPassCoefficients();
				break;
		}
	}

	void LowPassCoefficients() {
		float w0 = 2.0f * M_PI * _params.frequency / _fs;
		float cosw0 = cos(w0);
		float sinw0 = sin(w0);
		float alpha = sinw0/(2.0f*_params.q);

		_b0 = (1.0f - cosw0) / 2.0f;
		_b1 = 1.0f - cosw0;
		_b2 = (1.0f - cosw0) / 2.0f;
		_a0 = 1 + alpha;
		_a1 = -2.0f * cosw0;
		_a2 = 1 - alpha;
	}

	void HighPassCoefficients() {
		float w0 = 2.0f * M_PI * _params.frequency / _fs;
		float cosw0 = cos(w0);
		float sinw0 = sin(w0);
		float alpha = sinw0/(2.0f*_params.q);

		_b0 = (1.0f + cosw0) / 2.0f;
		_b1 = -(1.0f + cosw0);
		_b2 = (1.0f + cosw0) / 2.0f;
		_a0 = 1 + alpha;
		_a1 = -2.0f * cosw0;
		_a2 = 1 - alpha;
	}

	void BandPassCoefficients() {
		float w0 = 2.0f * M_PI * _params.frequency / _fs;
		float cosw0 = cos(w0);
		float sinw0 = sin(w0);
		float alpha = sinw0/(2.0f*_params.q);

		_b0 = alpha;
		_b1 = 0;
		_b2 = -alpha;
		_a0 = 1 + alpha;
		_a1 = -2.0f * cosw0;
		_a2 = 1 - alpha;
	}

	void NotchCoefficients() {
		float w0 = 2.0f * M_PI * _params.frequency / _fs;
		float cosw0 = cos(w0);
		float sinw0 = sin(w0);
		float alpha = sinw0/(2.0f*_params.q);

		_b0 = 1;
		_b1 = -2.0f * cosw0;
		_b2 = 1;
		_a0 = 1 + alpha;
		_a1 = -2.0f * cosw0;
		_a2 = 1 - alpha;
	}

	void AllPassCoefficients() {
		float w0 = 2.0f * M_PI * _params.frequency / _fs;
		float cosw0 = cos(w0);
		float sinw0 = sin(w0);
		float alpha = sinw0/(2.0f*_params.q);

		_b0 = 1 - alpha;
		_b1 = -2.0f * cosw0;
		_b2 = 1 + alpha;
		_a0 = 1 + alpha;
		_a1 = -2.0f * cosw0;
		_a2 = 1 - alpha;
	}
};

class DelaySoundChain : public SoundChainBase {
public:
	struct Parameters {
		float feedback = -0.5f;
		float delayTime = 0.15f;
	};

	DelaySoundChain() {}
	DelaySoundChain(Parameters &parameters) : _params(parameters ) {}
	~DelaySoundChain() {
		if (_buffer) delete[] _buffer;
	}

	Parameters GetParameters() {return _params;}
	void SetParameters(Parameters &parameters) {
		_params = parameters;

		_readIndex = _writeIndex - (_params.delayTime / BUFFER_DURATION * (float)_bufferSize);
		_readIndex = _readIndex % _bufferSize;
		_readIndex -= _readIndex % ReadSettings().Channels;
	}

private:
	const float BUFFER_DURATION = 2.0f;
	Parameters _params;
	float* _buffer = nullptr;
	int _bufferSize = 0;
	int _readIndex = 0;
	int _writeIndex = 0;

	void Reset() override {
		if (_buffer) delete[] _buffer;

		_bufferSize = ReadSettings().SampleRate * ReadSettings().Channels * BUFFER_DURATION;
		_buffer = new float[_bufferSize];
		for (int i = 0; i < _bufferSize; i++) {
			_buffer[i] = 0.0f;
		}

		_writeIndex = 0;
		_readIndex = _writeIndex - (_params.delayTime / BUFFER_DURATION * _bufferSize);
		_readIndex = _readIndex % _bufferSize;
		_readIndex -= _readIndex % ReadSettings().Channels;
	}

	void Process(float* buffPtr, int numberOfFrames) override {
		if (_buffer == nullptr) return;

		int numberOfChannels = ReadSettings().Channels;
		int numberOfSamples = numberOfFrames * numberOfChannels;
		for (int sample = 0; sample < numberOfSamples; sample++) {

			float value = _buffer[_readIndex];

			_buffer[_writeIndex] = buffPtr[sample] + _params.feedback * value;

			buffPtr[sample] += value;

			_readIndex++;
			_writeIndex++;
			if (_readIndex >= _bufferSize) _readIndex -= _bufferSize;
			if (_writeIndex >= _bufferSize) _writeIndex -= _bufferSize;
		}
	}
};
