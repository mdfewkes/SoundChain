#pragma once

#include "SoundChain.hpp"

// Wavetable sine osc
class PGSSinewaveSoundChain : public SoundChainBase {
public:
	PGSSinewaveSoundChain(double frequency) 
	: _freq(frequency) {}
	~PGSSinewaveSoundChain() {
		if (_waveTable) delete[] _waveTable;
	}

	void Reset() override {
		BuildTable();
	}

	double GetFrequency() { return _freq; }
	void SetFrequency(double frequency) {_freq = frequency; }

	float GetAmplitude() { return _amp; }
	void SetAmplitude(float amplitude) { _amp = amplitude; }

private:
	float* _waveTable = nullptr;
	int _tableSize = 0;
	double _positionInTable = 0;
	double _freq;
	float _amp = 0.25;

	void Process(float* buffPtr, int numberOfFrames) override {
		if (_waveTable == nullptr) return;

		for (int frame = 0; frame < numberOfFrames * _settings.Channels; frame += _settings.Channels) {
			float value = _amp * _waveTable[(int)_positionInTable];

			for (int channel = 0; channel < _settings.Channels; channel++) {
				buffPtr[frame + channel] += value;
			}

			_positionInTable += _freq;
			if (_positionInTable > _tableSize) _positionInTable -= _tableSize;
		}
	}

	void BuildTable() {
		if (_waveTable) delete[] _waveTable;

		_tableSize = _settings.SampleRate;

		_waveTable = new float[_tableSize];
		for (int position = 0; position < _tableSize; position++) {
			double phase = (double)position / (double)_tableSize;
			_waveTable[position] = sin(phase * 2 * M_PI);
		}
	}
};

// Sequel tests smaller table sizes
class PGSSinewaveSoundChain2 : public SoundChainBase {
public:
	PGSSinewaveSoundChain2(double frequency) 
	: _freq(frequency) {}
	~PGSSinewaveSoundChain2() {
		if (_waveTable != nullptr) delete[] _waveTable;
	}

	void Reset() override {
		BuildTable();
	}

	double GetFrequency() { return _freq; }
	void SetFrequency(double frequency) {_freq = frequency; }

	float GetAmplitude() { return _amp; }
	void SetAmplitude(float amplitude) { _amp = amplitude; }

private:
	float* _waveTable = nullptr;
	int _tableSize = 0;
	double _tableScale = 0.025;
	double _positionInTable = 0;
	double _freq;
	float _amp = 0.25;

	void Process(float* buffPtr, int numberOfFrames) override {
		if (_waveTable == nullptr) return;

		for (int frame = 0; frame < numberOfFrames * _settings.Channels; frame += _settings.Channels) {
			float value = _amp * _waveTable[(int)_positionInTable];

			for (int channel = 0; channel < _settings.Channels; channel++) {
				buffPtr[frame + channel] += value;
			}

			_positionInTable += _freq * _tableScale;
			if (_positionInTable > _tableSize) _positionInTable -= _tableSize;
		}
	}

	void BuildTable() {
		if (_waveTable) delete[] _waveTable;

		_tableScale = 2048.0 / _settings.SampleRate; // 2048 magic table target size

		_tableSize = _settings.SampleRate / (int)(1.0 / _tableScale);

		_waveTable = new float[_tableSize];
		for (int position = 0; position < _tableSize; position++) {
			double phase = (double)position / (double)_tableSize;
			_waveTable[position] = sin(phase * 2 * M_PI);
		}
	}
};

