#pragma once

#include "SoundChain.hpp"
#include <vector>


// Wavetable sine osc
class TBLAdditiveOscSoundChain : public SoundChainBase {
public:
	TBLAdditiveOscSoundChain(float frequency = 440.0f) : _baseFrequency(frequency) {}
	~TBLAdditiveOscSoundChain() {
		if (_waveTable) delete[] _waveTable;
	}

	float GetFrequency() {return _baseFrequency;}
	void SetFrequency(float frequency) {
		_baseFrequency = frequency;
		if (_oscFrequency.size() == 0) return;

		// float osc_mag = 0;
		numberOfActiveOsc = 0;
		for (int i = 0; i < NUMBER_OF_OSC; i++) {
			_oscFrequency.at(i) = _baseFrequency * _oscHarmonic.at(i);

			if (_oscFrequency.at(i) < ReadSettings().SampleRate / 2) {
				numberOfActiveOsc++;
				// osc_mag += _oscAmplitude.at(i);
			} else {
				break;
			}
		}
		// _amp = 1.0f / osc_mag * 0.5f;
	}

	void ImportTable(float* waveTable, size_t tableSize) {
		for (int i = NUMBER_OF_OSC-1; i >= 0; i--) {
			_oscPositionInTable.at(i) = (_oscPositionInTable.at(i) / (float)_tableSize) * (float)tableSize;
		}
		_waveTable = waveTable;
		_tableSize = tableSize;
		_tableScale = (float)_tableSize / (float)ReadSettings().SampleRate;
	}

private:
	const int DEFAULT_TABLE_SIZE = 64;
	const int NUMBER_OF_OSC = 32;
	int numberOfActiveOsc = 0;
	float* _waveTable = nullptr;
	int _tableSize = DEFAULT_TABLE_SIZE;
	float _tableScale;
	std::vector<float> _oscFrequency;
	std::vector<float> _oscAmplitude;
	std::vector<float> _oscPositionInTable;
	std::vector<float> _oscHarmonic;
	float _baseFrequency;
	float _amp = 0.5f;
	float _z1 = 0.0f;


	void Process(float* buffPtr, int numberOfFrames) override {
		if (_waveTable == nullptr) return;

		for (int frame = 0; frame < numberOfFrames * ReadSettings().Channels; frame += ReadSettings().Channels) {
			float value = 0.0f;

			for (int osc = 0; osc < numberOfActiveOsc; osc++) {
				value += _oscAmplitude.at(osc) * Lookup(_oscPositionInTable.at(osc));

				_oscPositionInTable.at(osc) += _oscFrequency.at(osc) * _tableScale;
				if (_oscPositionInTable.at(osc) >= _tableSize) _oscPositionInTable.at(osc) -= _tableSize;
			}

			value *= _amp;

			for (int channel = 0; channel < ReadSettings().Channels; channel++) {
				buffPtr[frame + channel] += value + _z1;
			}

			_z1 = value;
		}
	}

	void Reset() override {
		BuildTable();

		_oscFrequency.resize(NUMBER_OF_OSC, 0.0f);
		_oscAmplitude.resize(NUMBER_OF_OSC, 0.0f);
		_oscPositionInTable.resize(NUMBER_OF_OSC, 0.0f);
		_oscHarmonic.resize(NUMBER_OF_OSC, 0.0f);
		float osc_mag = 0;
		numberOfActiveOsc = 0;
		for (int i = 0; i < NUMBER_OF_OSC; i++) {
			_oscHarmonic.at(i) = i+1;
			_oscFrequency.at(i) = _baseFrequency * _oscHarmonic.at(i);
			_oscAmplitude.at(i) = 1.0f/(i+1);

			if (_oscFrequency.at(i) < ReadSettings().SampleRate / 2) {
				numberOfActiveOsc++;
				osc_mag += _oscAmplitude.at(i);
			}
		}
		_amp = 1.0f / osc_mag * 0.5f;
	}

	void BuildTable() {
		if (_waveTable) delete[] _waveTable;

		_tableScale = (float)_tableSize / (float)ReadSettings().SampleRate;

		_waveTable = new float[_tableSize];
		for (int position = 0; position < _tableSize; position++) {
			double phase = (double)position / (double)_tableSize;
			_waveTable[position] = sin(phase * 2 * M_PI);
		}
	}

	float Lookup(float position) {
		// return _waveTable[(int)position];

		float t = position - (int)position;
		float startPosition = (int)position;
		float endPosition = (int)position + 1;
		if (endPosition >= _tableSize) endPosition -= _tableSize;

		return (t * _waveTable[(int)endPosition]) + ((1.0f - t) * _waveTable[(int)startPosition]);
	}
};