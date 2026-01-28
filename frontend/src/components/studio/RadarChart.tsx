
import React from 'react';
import {
    Chart as ChartJS,
    RadialLinearScale,
    PointElement,
    LineElement,
    Filler,
    Tooltip,
    Legend,
} from 'chart.js';
import { Radar } from 'react-chartjs-2';

ChartJS.register(
    RadialLinearScale,
    PointElement,
    LineElement,
    Filler,
    Tooltip,
    Legend
);

interface RadarChartProps {
    baseline: Record<string, number>;
    proposal?: Record<string, number>;
}

export const RadarChart: React.FC<RadarChartProps> = ({ baseline, proposal }) => {
    const labels = ['Drug-likeness', 'Complexity', 'Stability', 'Solubility', 'Optimal Size'];

    // Helper to map record to array
    const mapScores = (scores: Record<string, number>) => [
        scores.drug_likeness || 0,
        scores.complexity || 0,
        scores.lipophilicity || scores.lipinski || 0,
        scores.solubility || 0,
        scores.size || 0,
    ];

    const datasets = [
        {
            label: 'Baseline',
            data: mapScores(baseline),
            backgroundColor: 'rgba(209, 213, 219, 0.2)',
            borderColor: 'rgb(156, 163, 175)',
            borderWidth: 1.5,
            pointBackgroundColor: 'rgb(156, 163, 175)',
            pointBorderColor: '#fff',
            pointRadius: 2,
        }
    ];

    if (proposal) {
        datasets.push({
            label: 'Proposal',
            data: mapScores(proposal),
            backgroundColor: 'rgba(59, 130, 246, 0.2)',
            borderColor: 'rgb(59, 130, 246)',
            borderWidth: 2,
            pointBackgroundColor: 'rgb(59, 130, 246)',
            pointBorderColor: '#fff',
            pointRadius: 3,
        });
    }

    const data = {
        labels: labels,
        datasets: datasets,
    };

    const options = {
        scales: {
            r: {
                angleLines: {
                    display: true,
                    color: 'rgba(0, 0, 0, 0.05)',
                },
                suggestedMin: 0,
                suggestedMax: 1,
                ticks: {
                    display: false,
                    stepSize: 0.2,
                },
                grid: {
                    color: 'rgba(0,0,0,0.05)'
                },
                pointLabels: {
                    font: {
                        size: 9,
                        weight: 'bold' as const,
                        family: 'Inter',
                    },
                    color: '#94a3b8',
                }
            },
        },
        plugins: {
            legend: {
                display: false,
            },
            tooltip: {
                enabled: true,
                backgroundColor: 'rgba(255, 255, 255, 0.9)',
                titleColor: '#1e293b',
                bodyColor: '#1e293b',
                borderColor: '#e2e8f0',
                borderWidth: 1,
                padding: 8,
                bodyFont: { size: 10 },
                displayColors: false
            }
        },
        maintainAspectRatio: false,
    };

    return (
        <div className="w-full h-40">
            <Radar data={data} options={options} />
        </div>
    );
};
