#pragma once

#include "SoundChain.hpp"
#include <atomic>
#include <cstdint>

struct SoundChainPlatformSettings {
	int SampleRate = 48000;
	int Channels = 2;
};

class SoundChainPlatform {
protected:

	virtual void Setup() {};
	virtual void Start() {};
	virtual void End() {};

public:
	SoundChainPlatform() {};
	~SoundChainPlatform() {};

	void Initialize(SoundChainPlatformSettings soundChainPlatformSettings) {
		if (_initialized) return;

		_settings = soundChainPlatformSettings;

		Setup();

		if (_previous) {
			_previous->Initialize(GetSoundChainSettings());
		}

		_samplesElapsed.store(0);

		Start();

		_initialized = true;
	};

	bool IsInitialized() {return _initialized;};

	void Terminate() {
		if (!_initialized) return;

		End();

		_initialized = false;
	};

	void FillBuffer(float* buffPtr, int numberOfFrames) {
		_previous->ReadSamples(buffPtr, numberOfFrames);
		_samplesElapsed.fetch_add(numberOfFrames);
	}

	SoundChainBase* GetPrevious() {return _previous;};
	void SetPrevious(SoundChainBase* previous) {_previous = previous;};

	SoundChainSettings GetSoundChainSettings() {
		SoundChainSettings settings = {_settings.SampleRate, _settings.Channels};
		return settings;
	};
	SoundChainPlatformSettings GetSoundChainPlatformSettings() {return _settings;};

	unsigned int GetSamplesElapsed() {return _samplesElapsed.load();};
	double GetTime() {return (double)_samplesElapsed.load() / (double)_settings.SampleRate;};

private:
	SoundChainPlatformSettings _settings;
	SoundChainBase* _previous = nullptr;
	bool _initialized = false;
	std::atomic<uint64_t> _samplesElapsed;
};