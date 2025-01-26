#pragma once

#include "SoundChain.hpp"

// TODO: Sequel with threaded buffer fill?
class DoubleBufferSoundChain : public SoundChainBase {
public:
	DoubleBufferSoundChain(size_t bufferFramesInSamples = 4096) 
	: _numberOfFrames(bufferFramesInSamples) {}
	~DoubleBufferSoundChain() {
		DestroyBuffers();
	}

	void Reset() override {
		DestroyBuffers();
		InitializeBuffers();

		if (_bufferInput) {
			_bufferInput->Initialize(_settings);
			FillSecondaryBuffer();
		}
	}

	SoundChainBase* GetBufferInput() {return _bufferInput;}
	void SetBufferInput(SoundChainBase* bufferInput) {
		_bufferInput = bufferInput;
	}

private:
	size_t _numberOfFrames;
	SoundChainBase* _bufferInput = nullptr;
	float* _bufferPointerPrimary = nullptr;
	float* _bufferPointerSecondary = nullptr;
	size_t _bufferSize = 0;
	size_t _bufferIndex = 0;

	void InitializeBuffers() {
		_bufferSize = _numberOfFrames * _settings.Channels;
		_bufferIndex = _bufferSize;

		_bufferPointerPrimary = new float[_bufferSize];
		_bufferPointerSecondary = new float[_bufferSize];
	}

	void DestroyBuffers() {
		if (_bufferPointerPrimary != nullptr) {
			delete[] _bufferPointerPrimary;
		}
		if (_bufferPointerSecondary != nullptr) {
			delete[] _bufferPointerSecondary;
		}
		_bufferSize = 0;
		_bufferIndex = 0;
	}

	void FillSecondaryBuffer() {
		if (_bufferInput == nullptr) return;
		
		_bufferInput->ReadSamples(_bufferPointerSecondary, _bufferSize/_settings.Channels);
	}

	void Process(float* buffPtr, int numberOfFrames) override {
		int numberOfSamples = numberOfFrames * _settings.Channels;
		for (int sample = 0; sample < numberOfSamples; sample++) {

			if (_bufferIndex >= _bufferSize) {
				// Swap buffers
				float* offhandBuffer = _bufferPointerPrimary;
				_bufferPointerPrimary = _bufferPointerSecondary;
				_bufferPointerSecondary = offhandBuffer;

				FillSecondaryBuffer();

				_bufferIndex = 0;
			}

			buffPtr[sample] += _bufferPointerPrimary[_bufferIndex];

			_bufferIndex++;
		}
	}
};

// TODO: Sequel with threaded buffer fill?
class MultiBufferSoundChain : public SoundChainBase {
public:
	MultiBufferSoundChain(size_t numberOfBuffers = 2, size_t bufferFramesInSamples = 4096)
	: _numberOfFrames(bufferFramesInSamples), _numberOfBuffers(numberOfBuffers) {}
	~MultiBufferSoundChain() {
		DestroyBuffers();
	}

	void Reset() override {
		DestroyBuffers();
		InitializeBuffers();

		if (_bufferInput) {
			_bufferInput->Initialize(_settings);
			FillLastBuffer();
		}
	}

	SoundChainBase* GetBufferInput() {return _bufferInput;}
	void SetBufferInput(SoundChainBase* bufferInput) {
		_bufferInput = bufferInput;
	}

private:
	size_t _numberOfFrames;
	SoundChainBase* _bufferInput = nullptr;
	float** _bufferPointers = nullptr;
	size_t _numberOfBuffers = 0;
	size_t _currentPointer = 0;
	size_t _bufferSize = 0;
	size_t _bufferIndex = 0;

	void InitializeBuffers() {
		_bufferSize = _numberOfFrames * _settings.Channels;
		_bufferIndex = _bufferSize;

		_bufferPointers = new float*[_numberOfBuffers];
		for (int i = 0; i < _numberOfBuffers; i ++) {
			_bufferPointers[i] = new float[_bufferSize];
		}
		_currentPointer = 0;
	}

	void DestroyBuffers() {
		if (_bufferPointers != nullptr) {
			for (int i = 0; i < _numberOfBuffers; i ++) {
				if (_bufferPointers[i] != nullptr) {
					delete[] _bufferPointers[i];
				}
			}

			delete[] _bufferPointers;
		}
	}

	void FillLastBuffer() {
		if (_bufferInput == nullptr) return;
		
		int bufferIndex = _currentPointer - 1;
		if (bufferIndex < 0) bufferIndex = _numberOfBuffers - 1;
		_bufferInput->ReadSamples(_bufferPointers[bufferIndex], _bufferSize/_settings.Channels);
	}

	void Process(float* buffPtr, int numberOfFrames) override {
		int numberOfSamples = numberOfFrames * _settings.Channels;
		for (int sample = 0; sample < numberOfSamples; sample++) {

			if (_bufferIndex >= _bufferSize) {
				_currentPointer++;
				if (_currentPointer >= _numberOfBuffers) {
					_currentPointer = 0;
				}

				FillLastBuffer();

				_bufferIndex = 0;
			}

			buffPtr[sample] += _bufferPointers[_currentPointer][_bufferIndex];

			_bufferIndex++;
		}
	}
};

// TODO: Circular buffer